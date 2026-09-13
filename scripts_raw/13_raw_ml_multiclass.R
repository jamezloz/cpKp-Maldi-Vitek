source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(data.table)
})

sample_group_split_multi <- function(y, groups, train_frac = 0.8, max_tries = 100) {
  ug <- unique(groups)
  for (i in seq_len(max_tries)) {
    tr_g <- sample(ug, size = floor(train_frac * length(ug)))
    tr <- groups %in% tr_g
    te <- !tr
    if (length(unique(y[tr])) == length(unique(y)) &&
        length(unique(y[te])) == length(unique(y))) {
      return(list(train = tr, test = te))
    }
  }
  NULL
}

balanced_acc <- function(y_true, y_pred) {
  mean(sapply(levels(y_true), function(cls) {
    idx <- y_true == cls
    if (!any(idx)) NA_real_ else mean(y_pred[idx] == cls)
  }), na.rm = TRUE)
}

macro_f1 <- function(y_true, y_pred) {
  classes <- levels(y_true)
  f1s <- sapply(classes, function(cls) {
    tp <- sum(y_true == cls & y_pred == cls)
    fp <- sum(y_true != cls & y_pred == cls)
    fn <- sum(y_true == cls & y_pred != cls)
    prec <- if (tp + fp == 0) 0 else tp / (tp + fp)
    rec  <- if (tp + fn == 0) 0 else tp / (tp + fn)
    if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
  })
  mean(f1s)
}

run_multiclass_cv <- function(X, y, groups, model_name = c("rf", "xgb"), n_iter = 100) {
  model_name <- match.arg(model_name)
  classes    <- levels(y)
  nclass     <- length(classes)
  out <- list()

  for (i in seq_len(n_iter)) {
    split <- sample_group_split_multi(y, groups)
    if (is.null(split)) next
    tr <- split$train; te <- split$test
    ytr <- y[tr]; yte <- y[te]

    if (model_name == "rf") {
      fit  <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500)
      pred <- predict(fit, X[te, , drop = FALSE])
    } else {
      ytr_num <- as.integer(ytr) - 1
      yte_num <- as.integer(yte) - 1
      dtr <- xgb.DMatrix(data = X[tr, , drop = FALSE], label = ytr_num)
      dte <- xgb.DMatrix(data = X[te, , drop = FALSE], label = yte_num)
      fit <- xgb.train(
        params = list(objective = "multi:softmax", num_class = nclass,
                      eval_metric = "mlogloss"),
        data = dtr, nrounds = 100, verbose = 0
      )
      pred_num <- predict(fit, dte)
      pred <- factor(classes[pred_num + 1], levels = classes)
    }

    out[[length(out) + 1]] <- data.table(
      iter     = i,
      model    = model_name,
      accuracy = mean(pred == yte),
      balanced_accuracy = balanced_acc(yte, pred),
      macro_f1 = macro_f1(yte, pred)
    )
  }
  rbindlist(out)
}

main <- function() {
  ds   <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_peak_strain_treatment_dataset.rds"))
  X    <- ds$matrix
  meta <- ds$metadata

  cat(sprintf("Features (sin prefiltro): %d\n", ncol(X)))

  out_dir <- file.path(RAW_MODULAR_OUT, "tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  schemes <- list(
    three_class = meta$mechanism %in% c("SEN", "KPC", "NDM"),
    four_class  = meta$mechanism %in% c("SEN", "KPC", "NDM", "KPC-NDM")
  )

  label_fn <- list(
    three_class = function(m) factor(m, levels = c("SEN", "KPC", "NDM")),
    four_class  = function(m) factor(m, levels = c("SEN", "KPC", "NDM", "KPC-NDM"))
  )

  all_results <- list()

  for (sc in names(schemes)) {
    keep <- schemes[[sc]]
    Xsc  <- X[keep, , drop = FALSE]
    y    <- label_fn[[sc]](meta$mechanism[keep])
    grps <- meta$strain_id[keep]
    cat("Scheme:", sc, "— n =", sum(keep), "— classes:", paste(levels(y), collapse="/"), "\n")

    for (mod in c("rf", "xgb")) {
      cat(" Model:", mod, "\n")
      res <- run_multiclass_cv(Xsc, y, grps, model_name = mod, n_iter = 100)
      res[, scheme := sc]
      all_results[[paste(sc, mod)]] <- res
    }
  }

  out <- rbindlist(all_results)
  fwrite(out, file.path(out_dir, "raw_multiclass_iterations.csv"))

  summary_dt <- out[, .(
    accuracy_mean  = mean(accuracy),
    accuracy_sd    = sd(accuracy),
    balanced_acc_mean = mean(balanced_accuracy),
    balanced_acc_sd   = sd(balanced_accuracy),
    macro_f1_mean  = mean(macro_f1),
    macro_f1_sd    = sd(macro_f1),
    n_iter         = .N
  ), by = .(scheme, model)]
  fwrite(summary_dt, file.path(out_dir, "raw_multiclass_summary.csv"))

  cat("\nSummary:\n")
  print(summary_dt)
  cat("Saved to:", out_dir, "\n")
}

set.seed(42)
main()
