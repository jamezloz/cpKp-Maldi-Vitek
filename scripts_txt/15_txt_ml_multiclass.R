source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(xgboost)
})

macro_balanced_accuracy <- function(truth, pred) {
  lv <- levels(truth)
  vals <- sapply(lv, function(cl) {
    idx <- truth == cl
    if (!any(idx)) return(NA_real_)
    mean(pred[idx] == cl)
  })
  mean(vals, na.rm = TRUE)
}

macro_f1 <- function(truth, pred) {
  lv <- levels(truth)
  vals <- sapply(lv, function(cl) {
    tp <- sum(truth == cl & pred == cl)
    fp <- sum(truth != cl & pred == cl)
    fn <- sum(truth == cl & pred != cl)

    if ((tp + fn) == 0 && (tp + fp) == 0) return(NA_real_)
    prec <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
  })
  mean(vals, na.rm = TRUE)
}

run_rf_multiclass <- function(X, y, groups, n_iter = 100) {
  ug <- unique(groups)
  rows <- list()
  for (i in seq_len(n_iter)) {
    tr_groups <- sample(ug, size = max(2, floor(0.8 * length(ug))))
    tr <- groups %in% tr_groups
    te <- !tr
    ytr <- droplevels(y[tr])
    yte <- droplevels(y[te])
    if (length(unique(ytr)) < length(levels(y))) next
    if (length(unique(yte)) < 2) next
    fit <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500)
    pred <- predict(fit, X[te, , drop = FALSE], type = "response")
    rows[[length(rows) + 1]] <- data.table(
      accuracy = mean(pred == yte),
      balanced_accuracy = macro_balanced_accuracy(yte, pred),
      macro_f1 = macro_f1(yte, pred)
    )
  }
  if (!length(rows)) return(data.table())
  rbindlist(rows)
}

run_xgb_multiclass <- function(X, y, groups, n_iter = 100) {
  ug <- unique(groups)
  rows <- list()
  y0 <- as.integer(y) - 1L
  nclass <- length(levels(y))
  for (i in seq_len(n_iter)) {
    tr_groups <- sample(ug, size = max(2, floor(0.8 * length(ug))))
    tr <- groups %in% tr_groups
    te <- !tr
    ytr <- droplevels(y[tr])
    yte <- y[te]
    if (length(unique(ytr)) < nclass) next
    if (length(unique(yte)) < 2) next
    ytr0 <- as.integer(ytr) - 1L
    dtrain <- xgb.DMatrix(data = as.matrix(X[tr, , drop = FALSE]), label = ytr0)
    dtest <- xgb.DMatrix(data = as.matrix(X[te, , drop = FALSE]))
    fit <- xgb.train(
      params = list(objective = "multi:softprob", eval_metric = "mlogloss", num_class = nclass),
      data = dtrain,
      nrounds = 100,
      verbose = 0
    )

    prob <- predict(fit, dtest)
    if (is.null(dim(prob))) prob <- matrix(prob, ncol = nclass, byrow = TRUE)
    stopifnot(ncol(prob) == nclass, nrow(prob) == sum(te))
    pred <- factor(levels(y)[max.col(prob, ties.method = "first")], levels = levels(y))
    rows[[length(rows) + 1]] <- data.table(
      accuracy = mean(pred == yte),
      balanced_accuracy = macro_balanced_accuracy(yte, pred),
      macro_f1 = macro_f1(yte, pred)
    )
  }
  if (!length(rows)) return(data.table())
  rbindlist(rows)
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X <- reps$intensity_raw

  scenarios <- list(
    three_class = meta$mechanism_raw %in% c("SEN", "KPC", "NDM", "KPC-OXA", "NDM-OXA"),
    four_class = meta$mechanism_raw %in% c("SEN", "KPC", "NDM", "KPC-NDM", "KPC-OXA", "NDM-OXA")
  )

  all_summary <- list()
  all_iter <- list()
  for (sc in names(scenarios)) {
    keep <- scenarios[[sc]]
    meta_sc <- copy(meta[keep])
    y <- if (sc == "three_class") {
      factor(fifelse(meta_sc$mechanism_raw == "SEN", "SEN", fifelse(meta_sc$mechanism_raw %in% c("KPC", "KPC-OXA"), "KPC", "NDM")), levels = c("SEN", "KPC", "NDM"))
    } else {
      factor(fifelse(meta_sc$mechanism_raw == "SEN", "SEN", fifelse(meta_sc$mechanism_raw %in% c("KPC", "KPC-OXA"), "KPC", fifelse(meta_sc$mechanism_raw %in% c("NDM", "NDM-OXA"), "NDM", "KPC-NDM"))), levels = c("SEN", "KPC", "NDM", "KPC-NDM"))
    }
    groups <- meta_sc$strain_id
    Xsc <- X[keep, , drop = FALSE]

    rf_res <- run_rf_multiclass(Xsc, y, groups, n_iter = 100)
    if (nrow(rf_res)) {
      rf_res[, `:=`(scenario = sc, model = "rf")]
      all_iter[[length(all_iter) + 1]] <- rf_res
      all_summary[[length(all_summary) + 1]] <- data.table(scenario = sc, model = "rf", n_classes = length(levels(y)), accuracy_mean = mean(rf_res$accuracy), balanced_accuracy_mean = mean(rf_res$balanced_accuracy), macro_f1_mean = mean(rf_res$macro_f1), n_iter_valid = nrow(rf_res))
    }

    xgb_res <- run_xgb_multiclass(Xsc, y, groups, n_iter = 100)
    if (nrow(xgb_res)) {
      xgb_res[, `:=`(scenario = sc, model = "xgb")]
      all_iter[[length(all_iter) + 1]] <- xgb_res
      all_summary[[length(all_summary) + 1]] <- data.table(scenario = sc, model = "xgb", n_classes = length(levels(y)), accuracy_mean = mean(xgb_res$accuracy), balanced_accuracy_mean = mean(xgb_res$balanced_accuracy), macro_f1_mean = mean(xgb_res$macro_f1), n_iter_valid = nrow(xgb_res))
    }
  }

  if (length(all_iter)) save_table(rbindlist(all_iter, fill = TRUE), "ml_multiclass_iterations")
  if (length(all_summary)) save_table(rbindlist(all_summary, fill = TRUE), "ml_multiclass_summary")

  cat("Saved multiclass ML outputs to:", TXT_MODULAR_OUT, "\n")
}

set.seed(42)
main()
