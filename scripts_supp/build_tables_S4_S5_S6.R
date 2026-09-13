source("scripts_txt/00_shared_utils.R")
suppressPackageStartupMessages(library(data.table))

OUT <- "resultados/PAPER_SUPP_TABLES"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dir_label <- function(cd, pos, neg) fifelse(cd > 0, pos, neg)

stab_res <- fread(file.path(TXT_MODULAR_TAB, "stability_selection_res_sen.csv"))
stab_kpc <- fread(file.path(TXT_MODULAR_TAB, "stability_selection_kpc_ndm.csv"))
shap_res <- fread(file.path(TXT_MODULAR_TAB, "shap_importance_res_sen.csv"))[order(-mean_abs_shap)]
shap_kpc <- fread(file.path(TXT_MODULAR_TAB, "shap_importance_kpc_ndm.csv"))[order(-mean_abs_shap)]
shap_res[, shap_rank := .I]; shap_kpc[, shap_rank := .I]

max_stab_in_range <- function(cmin, cmax, stab_dt) {
  vapply(seq_along(cmin), function(i) {
    v <- stab_dt[feature >= cmin[i] & feature <= cmax[i], selection_freq]
    if (!length(v)) NA_real_ else max(v)
  }, numeric(1))
}
shap_rank_in_range <- function(cmin, cmax, shap_dt) {
  vapply(seq_along(cmin), function(i) {
    v <- shap_dt[feature >= cmin[i] & feature <= cmax[i], shap_rank]
    if (!length(v)) NA_integer_ else min(v)
  }, integer(1))
}

build_s4_part <- function(collapsed_csv, stab_dt, shap_dt, pos_lab, neg_lab) {
  dt <- fread(collapsed_csv)
  dt <- dt[sig_wilcox == TRUE]
  dt[, stability := max_stab_in_range(cluster_min, cluster_max, stab_dt)]
  dt[, shap_rank := shap_rank_in_range(cluster_min, cluster_max, shap_dt)]
  dt[, direction := dir_label(cliff_delta, pos_lab, neg_lab)]
  setorder(dt, wilcox_fdr, wilcox_p)
  dt[, .(mz = cluster_center, cluster_min, cluster_max, wilcox_p, wilcox_fdr,
         fisher_p, fisher_fdr, auc_cv, cliff_delta, direction, stability, shap_rank)]
}

s4a <- build_s4_part(file.path(TXT_MODULAR_TAB, "final_biomarkers_res_sen_collapsed_10Da.csv"),
                     stab_res, shap_res, "CRKp", "CSKp")
s4b <- build_s4_part(file.path(TXT_MODULAR_TAB, "final_biomarkers_kpc_ndm_collapsed_10Da.csv"),
                     stab_kpc, shap_kpc, "KPC", "NDM")
fwrite(s4a, file.path(OUT, "S4_partA_crkp_cskp.csv"))
fwrite(s4b, file.path(OUT, "S4_partB_kpc_ndm.csv"))
cat(sprintf("S4 Part A (CRKp/CSKp): %d filas | Part B (KPC/NDM): %d filas\n", nrow(s4a), nrow(s4b)))

shap_res_raw <- fread("resultados/RAW_modular/tables/shap_importance_res_sen_raw.csv")[order(-mean_abs_shap)]
shap_kpc_raw <- fread("resultados/RAW_modular/tables/shap_importance_kpc_ndm_raw.csv")[order(-mean_abs_shap)]
shap_res_raw[, shap_rank := .I]; shap_kpc_raw[, shap_rank := .I]

nearest_shap_rank <- function(mz, shap_dt, tol = 2) {
  vapply(mz, function(m) {
    d <- abs(shap_dt$feature - m)
    i <- which.min(d)
    if (length(i) && d[i] <= tol) shap_dt$shap_rank[i] else NA_integer_
  }, integer(1))
}

build_s5_part <- function(tophat_csv, shap_dt, pos_lab, neg_lab) {
  dt <- fread(tophat_csv)
  dt <- dt[stability == 1]
  dt[, direction := dir_label(cliff_d, pos_lab, neg_lab)]
  dt[, shap_rank := nearest_shap_rank(mz, shap_dt)]
  setorder(dt, p_wilcox)
  dt[, .(mz, p_wilcox, fdr_wilcox, p_fisher, fdr_fisher, cliff_d, OR, pct_pos, pct_neg,
         sig_wilcox, sig_fisher, direction, shap_rank)]
}

s5a <- build_s5_part("resultados/RAW_modular/tables/tophat_biomarkers/tophat_strain_res_sen.csv",
                     shap_res_raw, "CRKp", "CSKp")
s5b <- build_s5_part("resultados/RAW_modular/tables/tophat_biomarkers/tophat_strain_kpc_ndm.csv",
                     shap_kpc_raw, "KPC", "NDM")
fwrite(s5a, file.path(OUT, "S5_partA_crkp_cskp.csv"))
fwrite(s5b, file.path(OUT, "S5_partB_kpc_ndm.csv"))
cat(sprintf("S5 Part A (CRKp/CSKp, stability=1.00): %d filas | Part B (KPC/NDM): %d filas\n", nrow(s5a), nrow(s5b)))

kpc <- fread("resultados/RAW_modular/tables/tophat_biomarkers/tophat_strain_kpc_ndm.csv")
s6 <- kpc[sig_fisher == TRUE & pct_neg < 0.05 & pct_pos > 0.4]
s6[, direction := "KPC (absent in NDM)"]
setorder(s6, -pct_pos)
s6 <- s6[, .(mz, p_wilcox, fdr_wilcox, p_fisher, fdr_fisher, OR, pct_pos, pct_neg, cliff_d, stability, direction)]
fwrite(s6, file.path(OUT, "S6_binary_kpc_specific.csv"))
cat(sprintf("S6 (marcadores binarios KPC-especificos, ausentes en NDM): %d filas\n", nrow(s6)))
print(s6)
