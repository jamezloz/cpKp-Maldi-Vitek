source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(randomForest)
  library(fastshap)
  library(ggplot2)
  library(data.table)
})

compute_shap_rf <- function(rf_obj, X_te, nsim = 100) {
  pfun <- function(object, newdata) predict(object, newdata, type = "prob")[, 2]
  shap <- fastshap::explain(rf_obj, X = X_te, pred_wrapper = pfun, nsim = nsim)
  as.matrix(shap)
}

run_shap_task <- function(X, y, task_label, pos_label, groups = NULL, n_iter = 10, nsim = 50) {
  set.seed(42)
  if (is.null(groups)) groups <- rownames(X)
  ug <- unique(groups)
  shap_list <- list()

  acc <- setNames(numeric(ncol(X)), colnames(X))
  cnt <- setNames(integer(ncol(X)), colnames(X))
  for (i in seq_len(n_iter)) {
    tr_groups <- sample(ug, size = floor(0.8 * length(ug)))
    tr <- groups %in% tr_groups
    te <- !tr
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) next

    ytr <- factor(y[tr])
    pv  <- apply(X[tr, , drop = FALSE], 2,
                 function(v) tryCatch(wilcox.test(v ~ ytr)$p.value, error = function(e) NA_real_))
    keep <- which(!is.na(pv) & p.adjust(pv, "fdr") < 0.05)
    if (length(keep) < 5) keep <- order(pv, na.last = TRUE)[seq_len(min(20, ncol(X)))]

    fit <- randomForest(x = X[tr, keep, drop = FALSE], y = ytr, ntree = 500)
    shap_mat <- compute_shap_rf(fit, X[te, keep, drop = FALSE], nsim = nsim)
    mabs <- colMeans(abs(shap_mat))
    acc[colnames(X)[keep]] <- acc[colnames(X)[keep]] + mabs
    cnt[colnames(X)[keep]] <- cnt[colnames(X)[keep]] + 1L
  }

  sel <- cnt > 0
  importance <- data.table(
    feature       = names(acc)[sel],
    mean_abs_shap = acc[sel] / cnt[sel],
    n_iter_sel    = cnt[sel],
    task          = task_label
  )[order(-mean_abs_shap)]
  importance
}

main <- function() {
  ds  <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_peak_strain_treatment_dataset.rds"))
  X   <- ds$matrix
  meta <- ds$metadata
  rownames(X) <- meta$sample

  cat(sprintf("  Features (sin prefiltro): %d\n", ncol(X)))

  out_dir <- file.path(RAW_MODULAR_OUT, "tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()

  keep_rs <- meta$label_res_sen %in% c("RES", "SEN")
  cat("Running SHAP: RES vs SEN\n")
  imp_rs <- run_shap_task(X[keep_rs, , drop = FALSE],
                          meta$label_res_sen[keep_rs], "res_sen", "RES",
                          groups = meta$strain_id[keep_rs])
  results[["res_sen"]] <- imp_rs
  fwrite(imp_rs, file.path(out_dir, "shap_importance_res_sen_raw.csv"))

  keep_kn <- meta$mechanism %in% c("KPC", "NDM")
  cat("Running SHAP: KPC vs NDM\n")
  imp_kn <- run_shap_task(X[keep_kn, , drop = FALSE],
                          meta$mechanism[keep_kn], "kpc_ndm", "KPC",
                          groups = meta$strain_id[keep_kn])
  results[["kpc_ndm"]] <- imp_kn
  fwrite(imp_kn, file.path(out_dir, "shap_importance_kpc_ndm_raw.csv"))

  for (task in names(results)) {
    top20 <- results[[task]][1:min(20, .N)]
    p <- ggplot(top20, aes(reorder(feature, mean_abs_shap), mean_abs_shap)) +
      geom_col(fill = "#4DBBD5") +
      coord_flip() +
      theme_bw() +
      labs(title = paste("SHAP importance (RAW strain-mean) —", task),
           x = "m/z feature", y = "Mean |SHAP|")
    ggsave(file.path(RAW_MODULAR_OUT, "figures",
                     paste0("shap_importance_", task, "_raw.png")),
           p, width = 7, height = 6)
  }

  cat("Saved RAW SHAP outputs to:", out_dir, "\n")
}

main()
