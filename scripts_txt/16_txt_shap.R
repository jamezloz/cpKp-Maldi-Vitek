suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(fastshap)
  library(pROC)
})

TXT_MODULAR_OBJ <- "resultados/TXT_modular/objects"
TXT_MODULAR_TAB <- "resultados/TXT_modular/tables"
TXT_MODULAR_MOD <- "resultados/TXT_modular/models"
dir.create(TXT_MODULAR_MOD, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

strain_split <- function(meta, y_col, train_frac = 0.8) {
  strains  <- unique(meta$strain_id)
  tr_str   <- sample(strains, floor(train_frac * length(strains)))
  list(train = meta$strain_id %in% tr_str,
       test  = !meta$strain_id %in% tr_str)
}

train_rf <- function(X_tr, y_tr, ntree = 500) {

  pvals <- apply(X_tr, 2, function(v) tryCatch(wilcox.test(v ~ y_tr)$p.value, error = function(e) NA_real_))
  adj_p <- p.adjust(pvals, "fdr")
  keep  <- !is.na(adj_p) & adj_p < 0.05
  if (sum(keep) < 5) keep <- order(pvals, na.last = TRUE)[seq_len(min(20, ncol(X_tr)))]
  X_sel <- X_tr[, keep, drop = FALSE]
  rf    <- randomForest(x = X_sel, y = y_tr, ntree = ntree, importance = TRUE)
  list(model = rf, features = colnames(X_sel))
}

compute_shap <- function(rf_obj, X_te, y_te, nsim = 100) {
  pfun <- function(object, newdata)
    predict(object, newdata, type = "prob")[, levels(y_te)[2]]
  shap <- fastshap::explain(rf_obj, X = X_te, pred_wrapper = pfun, nsim = nsim,
                            adjust = TRUE)
  as.matrix(shap)
}

run_task <- function(label, X, meta, y, ntree = 500, nsim = 100) {
  cat(sprintf("\n── %s ──\n", label))

  sp   <- strain_split(meta, train_frac = 0.8)
  X_tr <- X[sp$train, , drop = FALSE];  y_tr <- y[sp$train]
  X_te <- X[sp$test,  , drop = FALSE];  y_te <- y[sp$test]

  cat(sprintf("  Train: %d samples | Test: %d samples\n", nrow(X_tr), nrow(X_te)))

  obj <- train_rf(X_tr, y_tr, ntree = ntree)
  rf  <- obj$model
  feats <- obj$features

  probs <- predict(rf, X_te[, feats, drop = FALSE], type = "prob")[, levels(y)[2]]
  auc   <- as.numeric(pROC::auc(pROC::roc(as.numeric(y_te == levels(y)[2]), probs, quiet = TRUE)))
  cat(sprintf("  Test AUC: %.3f  |  Features selected: %d\n", auc, length(feats)))

  cat(sprintf("  Computing SHAP (nsim=%d)...\n", nsim))
  X_te_sel <- X_te[, feats, drop = FALSE]
  shap_mat  <- compute_shap(rf, X_te_sel, y_te, nsim = nsim)

  mean_abs  <- colMeans(abs(shap_mat))
  imp_dt    <- data.table(
    feature      = names(mean_abs),
    mean_abs_shap = mean_abs,
    rf_importance = importance(rf, type = 1)[names(mean_abs), "MeanDecreaseAccuracy"]
  )[order(-mean_abs_shap)]

  sign_vec <- sapply(colnames(shap_mat), function(f) {
    cor(X_te_sel[, f], shap_mat[, f], method = "spearman")
  })
  imp_dt[, direction := fifelse(sign_vec[feature] > 0,
                                levels(y)[2], levels(y)[1])]

  shap_long <- as.data.table(shap_mat)
  shap_long[, sample_id := seq_len(nrow(shap_mat))]
  shap_long[, true_label := as.character(y_te)]
  shap_long[, pred_prob  := probs]

  list(imp = imp_dt, shap_matrix = shap_mat, X_te_sel = X_te_sel,
       shap_long = shap_long, auc = auc, model = rf, features = feats)
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X    <- reps$intensity_raw

  keep_res <- !is.na(meta$label_res_sen)
  y_res    <- factor(meta$label_res_sen[keep_res], levels = c("SEN", "RES"))
  res <- run_task("RES vs SEN", X[keep_res, ], meta[keep_res], y_res)

  fwrite(res$imp, file.path(TXT_MODULAR_TAB, "shap_importance_res_sen.csv"))
  shap_out      <- as.data.table(res$shap_matrix)
  shap_out[, true_label := res$shap_long$true_label]
  shap_out[, pred_prob  := res$shap_long$pred_prob]
  fval_out      <- as.data.table(res$X_te_sel)
  fwrite(shap_out,  file.path(TXT_MODULAR_TAB, "shap_matrix_res_sen.csv"))
  fwrite(fval_out,  file.path(TXT_MODULAR_TAB, "shap_feature_values_res_sen.csv"))
  saveRDS(res$model, file.path(TXT_MODULAR_MOD, "rf_res_sen.rds"))
  cat(sprintf("  RES/SEN: top 10 SHAP features:\n"))
  print(res$imp[1:10, .(feature, mean_abs_shap, direction)])

  keep_mech <- meta$mechanism %in% c("KPC", "NDM")
  y_mech    <- factor(meta$mechanism[keep_mech], levels = c("NDM", "KPC"))
  kn <- run_task("KPC vs NDM", X[keep_mech, ], meta[keep_mech], y_mech)

  fwrite(kn$imp, file.path(TXT_MODULAR_TAB, "shap_importance_kpc_ndm.csv"))
  shap_out2 <- as.data.table(kn$shap_matrix)
  shap_out2[, true_label := kn$shap_long$true_label]
  shap_out2[, pred_prob  := kn$shap_long$pred_prob]
  fval_out2 <- as.data.table(kn$X_te_sel)
  fwrite(shap_out2, file.path(TXT_MODULAR_TAB, "shap_matrix_kpc_ndm.csv"))
  fwrite(fval_out2, file.path(TXT_MODULAR_TAB, "shap_feature_values_kpc_ndm.csv"))
  saveRDS(kn$model, file.path(TXT_MODULAR_MOD, "rf_kpc_ndm.rds"))
  cat(sprintf("  KPC/NDM: top 10 SHAP features:\n"))
  print(kn$imp[1:10, .(feature, mean_abs_shap, direction)])

  cat("\nSaved SHAP outputs to:", TXT_MODULAR_TAB, "\n")
}

main()
