source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(pROC)
  library(data.table)
})

collapse_lab_replicates <- function(mat, meta) {
  key <- paste(meta$strain_id, meta$treatment, sep = "__")
  idx <- split(seq_len(nrow(meta)), key)
  m <- t(vapply(idx, function(i) colMeans(mat[i, , drop = FALSE]), numeric(ncol(mat))))
  md <- rbindlist(lapply(idx, function(i) {
    g <- meta[i, ]
    data.table(strain_id = g$strain_id[1], treatment = g$treatment[1],
               mechanism = g$mechanism[1], label_res_sen = g$label_res_sen[1])
  }))
  list(matrix = m, metadata = md)
}

univariate_auc_cv <- function(mat, y, n_folds = 5) {
  folds <- sample(rep(seq_len(n_folds), length.out = length(y)))
  rows <- lapply(seq_len(ncol(mat)), function(j) {
    probs <- rep(NA_real_, length(y))
    for (k in seq_len(n_folds)) {
      tr <- folds != k
      te <- folds == k
      df_tr <- data.frame(y = y[tr], x = mat[tr, j])
      df_te <- data.frame(x = mat[te, j])
      fit <- glm(y ~ x, data = df_tr, family = binomial())
      probs[te] <- predict(fit, newdata = df_te, type = "response")
    }
    auc_val <- as.numeric(pROC::auc(pROC::roc(y, probs, quiet = TRUE)))
    data.table(feature = colnames(mat)[j], auc_cv = auc_val)
  })
  rbindlist(rows)[order(-auc_cv)]
}

cliff_delta_vec <- function(mat, y) {
  lv <- levels(y)
  rows <- lapply(seq_len(ncol(mat)), function(j) {
    a <- mat[y == lv[2], j]; b <- mat[y == lv[1], j]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    n1 <- length(a); n2 <- length(b)
    cd <- if (!n1 || !n2) NA_real_ else {
      r <- rank(c(a, b)); 2 * ((sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)) - 1
    }
    data.table(feature = colnames(mat)[j], cliff_delta = cd)
  })
  rbindlist(rows)
}

pairwise_wilcox_simple <- function(mat, y) {
  rows <- lapply(seq_len(ncol(mat)), function(j) {
    p <- tryCatch(wilcox.test(mat[, j] ~ y)$p.value, error = function(e) NA_real_)
    data.table(feature = colnames(mat)[j], p_value = p)
  })
  out <- rbindlist(rows)
  out[, fdr := p.adjust(p_value, method = "fdr")]
  out[order(fdr, p_value)]
}

pairwise_fisher_simple <- function(mat_bin, y) {
  rows <- lapply(seq_len(ncol(mat_bin)), function(j) {
    tbl <- table(mat_bin[, j] > 0, y)
    if (nrow(tbl) < 2 || ncol(tbl) < 2) p <- NA_real_ else p <- fisher.test(tbl)$p.value
    data.table(feature = colnames(mat_bin)[j], p_value = p)
  })
  out <- rbindlist(rows)
  out[, fdr := p.adjust(p_value, method = "fdr")]
  out[order(fdr, p_value)]
}

assemble_task <- function(wil, auc_dt, fisher_dt, cliff_dt, taskname) {
  dt <- merge(wil[, .(feature, wilcox_p = p_value, wilcox_fdr = fdr)],
              fisher_dt[, .(feature, fisher_p = p_value, fisher_fdr = fdr)],
              by = "feature", all = TRUE)
  dt <- merge(dt, auc_dt, by = "feature", all = TRUE)
  dt <- merge(dt, cliff_dt, by = "feature", all = TRUE)
  dt <- dt[!(is.na(wilcox_p) & is.na(fisher_p))]
  dt[, sig_wilcox := !is.na(wilcox_fdr) & wilcox_fdr < 0.05]
  dt[, sig_fisher := !is.na(fisher_fdr) & fisher_fdr < 0.05]
  dt[, sig_ambos  := sig_wilcox & sig_fisher]
  dt[, tarea := taskname]
  setorder(dt, -sig_ambos, wilcox_fdr, na.last = TRUE)
  dt[]
}

main <- function() {
  set.seed(42)
  meta <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))

  binR    <- reps$intensity_raw > 0
  res_all <- meta$label_res_sen == "RES"
  prev_max <- pmax(colMeans(binR[res_all, , drop = FALSE]),
                   colMeans(binR[!res_all, , drop = FALSE]))
  keep_prev <- prev_max >= 0.05
  ds  <- collapse_lab_replicates(reps$intensity_raw[, keep_prev, drop = FALSE], meta)
  dsb <- collapse_lab_replicates(((reps$binary > 0) * 1)[, keep_prev, drop = FALSE], meta)
  X   <- ds$matrix; Xb <- dsb$matrix; md <- ds$metadata
  cat(sprintf("  Features canonicas (prevalencia maxima por clase >= 5%%): %d\n", ncol(X)))

  y_res <- factor(md$label_res_sen, levels = c("SEN", "RES"))
  final_res <- assemble_task(pairwise_wilcox_simple(X, y_res),
                             univariate_auc_cv(X, y_res, n_folds = 5),
                             pairwise_fisher_simple(Xb, y_res),
                             cliff_delta_vec(X, y_res), "res_sen")

  keep_mech <- md$mechanism %in% c("KPC", "NDM")
  y_mech <- factor(md$mechanism[keep_mech], levels = c("NDM", "KPC"))
  final_mech <- assemble_task(pairwise_wilcox_simple(X[keep_mech, , drop = FALSE], y_mech),
                              univariate_auc_cv(X[keep_mech, , drop = FALSE], y_mech, n_folds = 5),
                              pairwise_fisher_simple(Xb[keep_mech, , drop = FALSE], y_mech),
                              cliff_delta_vec(X[keep_mech, , drop = FALSE], y_mech), "kpc_ndm")

  sig_rs <- final_res[sig_wilcox | sig_fisher, feature]
  sig_kn <- final_mech[sig_wilcox | sig_fisher, feature]
  final_res[, tambien_sig_KPC_NDM := feature %in% sig_kn]
  final_mech[, tambien_sig_RES_SEN := feature %in% sig_rs]

  final_res_collapsed  <- collapse_ranked_biomarkers(final_res,  feature_col = "feature", score_cols = c("wilcox_fdr", "wilcox_p"), window_da = 10)
  final_mech_collapsed <- collapse_ranked_biomarkers(final_mech, feature_col = "feature", score_cols = c("wilcox_fdr", "wilcox_p"), window_da = 10)

  save_table(final_res, "final_biomarkers_res_sen_uncollapsed")
  save_table(final_res_collapsed, "final_biomarkers_res_sen_collapsed_10Da")
  save_table(final_mech, "final_biomarkers_kpc_ndm_uncollapsed")
  save_table(final_mech_collapsed, "final_biomarkers_kpc_ndm_collapsed_10Da")

  cat(sprintf("\nRES/SEN: %d testeables | sig_wilcox %d | sig_fisher %d | sig_ambos %d\n",
      nrow(final_res), final_res[, sum(sig_wilcox)], final_res[, sum(sig_fisher)], final_res[, sum(sig_ambos)]))
  cat(sprintf("KPC/NDM: %d testeables | sig_wilcox %d | sig_fisher %d | sig_ambos %d\n",
      nrow(final_mech), final_mech[, sum(sig_wilcox)], final_mech[, sum(sig_fisher)], final_mech[, sum(sig_ambos)]))
  cat(sprintf("Biomarcadores de resistencia que ademas separan mecanismo: %d\n",
      final_res[sig_ambos & tambien_sig_KPC_NDM, .N]))
  cat("Saved final biomarker tables to:", TXT_MODULAR_OUT, "\n")
}

main()
