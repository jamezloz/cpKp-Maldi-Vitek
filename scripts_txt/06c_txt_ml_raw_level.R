source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(randomForest)
  library(xgboost)
  library(pROC)
})

cv_ml_raw <- function(X, y, groups, model_name = c("rf", "xgb"), n_iter = 100) {
  model_name <- match.arg(model_name)
  ug        <- unique(groups)
  pos_label <- levels(y)[2]
  neg_label <- levels(y)[1]
  out <- list()

  for (i in seq_len(n_iter)) {
    tr_groups <- sample(ug, size = floor(0.8 * length(ug)))
    tr <- groups %in% tr_groups
    te <- !tr
    ytr <- y[tr]; yte <- y[te]
    if (length(unique(ytr)) < 2 || length(unique(yte)) < 2) next

    if (model_name == "rf") {
      fit   <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500)
      probs <- predict(fit, X[te, , drop = FALSE], type = "prob")[, pos_label]
      pred  <- predict(fit, X[te, , drop = FALSE], type = "response")
    } else {
      dtr  <- xgb.DMatrix(data = X[tr, , drop = FALSE], label = as.numeric(ytr == pos_label))
      dte  <- xgb.DMatrix(data = X[te, , drop = FALSE], label = as.numeric(yte == pos_label))
      fit  <- xgb.train(
        params = list(objective = "binary:logistic", eval_metric = "auc"),
        data = dtr, nrounds = 100, verbose = 0
      )
      probs <- predict(fit, dte)
      pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = levels(y))
    }

    tab  <- table(yte, pred)
    sens <- tab[pos_label, pos_label] / sum(tab[pos_label, ])
    spec <- tab[neg_label, neg_label] / sum(tab[neg_label, ])
    auc  <- as.numeric(pROC::auc(pROC::roc(as.numeric(yte == pos_label), probs, quiet = TRUE)))

    out[[length(out) + 1]] <- data.table(
      iter              = i,
      model             = model_name,
      level             = "raw_spectrum",
      accuracy          = mean(pred == yte),
      balanced_accuracy = mean(c(sens, spec), na.rm = TRUE),
      auc               = auc,
      n_train           = sum(tr),
      n_test            = sum(te)
    )
  }
  rbindlist(out)
}

main <- function() {
  meta_raw <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds"))
  mat_raw  <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))

  bin_raw     <- mat_raw > 0
  res_i       <- meta_raw$label_res_sen == "RES"
  prev_max    <- pmax(colMeans(bin_raw[res_i, , drop = FALSE]),
                      colMeans(bin_raw[!res_i, , drop = FALSE]))
  keep_feat   <- prev_max >= 0.05
  X <- mat_raw[, keep_feat, drop = FALSE]
  cat("Features after filtering:", ncol(X), "\n")

  tasks <- list(
    res_sen = list(
      keep = meta_raw$label_res_sen %in% c("RES", "SEN"),
      y_fn = function(m) factor(m$label_res_sen, levels = c("SEN", "RES"))
    ),
    kpc_ndm = list(
      keep = meta_raw$mechanism %in% c("KPC", "NDM"),
      y_fn = function(m) factor(m$mechanism, levels = c("NDM", "KPC"))
    )
  )

  all_res <- list()

  for (task in names(tasks)) {
    keep <- tasks[[task]]$keep
    Xk   <- X[keep, , drop = FALSE]
    y    <- tasks[[task]]$y_fn(meta_raw[keep])
    grps <- meta_raw$strain_id[keep]
    cat("Task:", task, "| n spectra:", sum(keep),
        "| n strains:", length(unique(grps)), "\n")

    for (mod in c("rf", "xgb")) {
      cat(" model:", mod, "\n")
      res <- cv_ml_raw(Xk, y, grps, model_name = mod, n_iter = 100)
      res[, task := task]
      all_res[[paste(task, mod)]] <- res
    }
  }

  out <- rbindlist(all_res)
  save_table(out, "ml_raw_level_iterations")

  summary_dt <- out[, .(
    balanced_accuracy_mean = mean(balanced_accuracy),
    balanced_accuracy_sd   = sd(balanced_accuracy),
    auc_mean               = mean(auc),
    auc_sd                 = sd(auc),
    n_iter                 = .N
  ), by = .(task, model, level)]
  save_table(summary_dt, "ml_raw_level_summary")

  combos <- expand.grid(task = c("res_sen", "kpc_ndm"), model = c("rf", "xgb"),
                        stringsAsFactors = FALSE)
  agg <- rbindlist(lapply(seq_len(nrow(combos)), function(k) {
    tk <- combos$task[k]; md <- combos$model[k]
    f  <- file.path(TXT_MODULAR_TAB,
                    sprintf("ml_internal_%s_%s_iterations.csv", tk, md))
    if (!file.exists(f)) {
      warning("FALTA el agregado ", basename(f), ": ", tk, "/", md,
              " no entra en la comparacion", call. = FALSE, immediate. = TRUE)
      return(NULL)
    }
    d <- fread(f)
    d[, `:=`(level = "strain_protocolo", task = tk, model = md)]
    d[, .(task, model, level, balanced_accuracy, auc)]
  }))

  comp <- rbindlist(list(
    out[, .(task, model, level, balanced_accuracy, auc)],
    agg
  ), fill = TRUE)

  comp_summary <- comp[, .(
    balanced_accuracy_mean = mean(balanced_accuracy),
    balanced_accuracy_sd   = sd(balanced_accuracy),
    auc_mean               = mean(auc),
    auc_sd                 = sd(auc)
  ), by = .(task, model, level)]
  save_table(comp_summary, "ml_raw_vs_aggregated_comparison")

  cat("\n=== Raw spectrum vs strain-mean comparison ===\n")
  print(comp_summary[order(task, model, level)])
  cat("Saved to:", TXT_MODULAR_TAB, "\n")
}

set.seed(42)
main()
