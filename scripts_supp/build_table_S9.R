suppressPackageStartupMessages(library(data.table))
OUT <- "resultados/PAPER_SUPP_TABLES"
TXT_TAB <- "resultados/TXT_modular/tables"
RAW_TAB <- "resultados/RAW_modular/tables"

build_txt <- function(shap_csv, stab_csv, pos_lab, neg_lab, n_top = 30) {
  s <- fread(shap_csv)[order(-mean_abs_shap)][seq_len(n_top)]
  st <- fread(stab_csv)
  s <- merge(s, st, by = "feature", all.x = TRUE)
  setorder(s, -mean_abs_shap)
  s[, rank := .I]
  s[, stable := fifelse(!is.na(selection_freq) & selection_freq >= 1.0, "Yes", "---")]
  s[, direction := fifelse(direction == pos_lab, pos_lab, neg_lab)]
  s[, .(rank, mz = feature, mean_abs_shap = round(mean_abs_shap, 5), direction, stable)]
}

txt_rs  <- build_txt(file.path(TXT_TAB, "shap_importance_res_sen.csv"),
                     file.path(TXT_TAB, "stability_selection_res_sen.csv"), "RES", "SEN")
txt_rs[, direction := fifelse(direction == "RES", "CRKp", "CSKp")]
txt_kn  <- build_txt(file.path(TXT_TAB, "shap_importance_kpc_ndm.csv"),
                     file.path(TXT_TAB, "stability_selection_kpc_ndm.csv"), "KPC", "NDM")

fwrite(txt_rs, file.path(OUT, "S9_txt_res_sen.csv"))
fwrite(txt_kn, file.path(OUT, "S9_txt_kpc_ndm.csv"))

build_raw <- function(shap_csv, tophat_csv, pos_lab, neg_lab, n_top = 30, tol = 1.5) {
  s <- fread(shap_csv)[order(-mean_abs_shap)][seq_len(n_top)]
  th <- fread(tophat_csv)
  s[, rank := .I]
  s[, `:=`(cliff_d = NA_real_, stability = NA_real_)]
  for (i in seq_len(nrow(s))) {
    d <- abs(th$mz - s$feature[i])
    j <- which.min(d)
    if (length(j) && d[j] <= tol) {
      s$cliff_d[i]   <- th$cliff_d[j]
      s$stability[i] <- th$stability[j]
    }
  }
  s[, direction := fifelse(is.na(cliff_d), "---", fifelse(cliff_d > 0, pos_lab, neg_lab))]
  s[, stable := fifelse(!is.na(stability) & stability >= 1.0, "Yes", "---")]
  s[, .(rank, mz = round(feature, 1), mean_abs_shap = round(mean_abs_shap, 5), direction, stable)]
}

raw_rs <- build_raw(file.path(RAW_TAB, "shap_importance_res_sen_raw.csv"),
                    file.path(RAW_TAB, "tophat_biomarkers/tophat_strain_res_sen.csv"), "CRKp", "CSKp")
raw_kn <- build_raw(file.path(RAW_TAB, "shap_importance_kpc_ndm_raw.csv"),
                    file.path(RAW_TAB, "tophat_biomarkers/tophat_strain_kpc_ndm.csv"), "KPC", "NDM")

fwrite(raw_rs, file.path(OUT, "S9_raw_res_sen.csv"))
fwrite(raw_kn, file.path(OUT, "S9_raw_kpc_ndm.csv"))

cat("TXT res_sen:", nrow(txt_rs), "| TXT kpc_ndm:", nrow(txt_kn),
    "| RAW res_sen:", nrow(raw_rs), "| RAW kpc_ndm:", nrow(raw_kn), "\n")
print(head(txt_rs, 8))
print(head(raw_kn, 8))
