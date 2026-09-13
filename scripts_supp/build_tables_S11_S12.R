suppressPackageStartupMessages(library(data.table))

OUT <- "resultados/PAPER_SUPP_TABLES"
PC  <- "resultados/RAW_modular/tables/panel_complementarity"
TXT <- "resultados/TXT_modular/tables"

kpc_steps <- fread(file.path(PC, "kpc_greedy_steps.csv"))
kpc_steps[, panel := paste0("+", mz, " Da")]
kpc_steps[paso == 0, panel := paste0(mz, " Da (seed)")]
fwrite(kpc_steps[, .(step = paso, mz, panel, sensitivity = sens, specificity = spec, new_TP = nuevos)],
       file.path(OUT, "S11_kpc_greedy_steps.csv"))

ndm_steps <- fread(file.path(PC, "ndm_greedy_steps.csv"))
ndm_steps[, panel := paste0("+", mz, " Da")]
fwrite(ndm_steps[, .(step = paso, mz, panel, sensitivity = sens, specificity = spec, new_TP = nuevos)],
       file.path(OUT, "S11_ndm_greedy_steps.csv"))

ndm_single <- fread(file.path(PC, "ndm_single_markers.csv"))[order(-youden)][1:5]
fwrite(ndm_single, file.path(OUT, "S11_ndm_top5_single.csv"))

oob_kpc <- fread(file.path(PC, "bootstrap_oob_KPC.csv"))
oob_ndm <- fread(file.path(PC, "bootstrap_oob_NDM.csv"))
oob_summary <- data.table(
  target = c("KPC", "NDM"),
  n_boot = c(nrow(oob_kpc), nrow(oob_ndm)),
  median_k = c(median(oob_kpc$k), median(oob_ndm$k)),
  k_range = c(paste(range(oob_kpc$k), collapse="-"), paste(range(oob_ndm$k), collapse="-")),
  median_sens = c(median(oob_kpc$sens), median(oob_ndm$sens)),
  sens_lo = c(quantile(oob_kpc$sens, .025), quantile(oob_ndm$sens, .025)),
  sens_hi = c(quantile(oob_kpc$sens, .975), quantile(oob_ndm$sens, .975)),
  median_spec = c(median(oob_kpc$spec), median(oob_ndm$spec)),
  spec_lo = c(quantile(oob_kpc$spec, .025), quantile(oob_ndm$spec, .025)),
  spec_hi = c(quantile(oob_kpc$spec, .975), quantile(oob_ndm$spec, .975))
)
fwrite(oob_summary, file.path(OUT, "S11_oob_summary.csv"))
cat("S11: KPC steps", nrow(kpc_steps), "| NDM steps", nrow(ndm_steps),
    "| NDM top5 single", nrow(ndm_single), "| OOB summary 2 rows\n")
print(oob_summary)

pairs <- fread(file.path(TXT, "kpc_ndm_platform_concordance_pairs.csv"))
setorder(pairs, strain_id, treatment)
fwrite(pairs[, .(strain_id, treatment, truth, call_site_A = call_A, call_site_M = call_M,
                 concordant = concordante)],
       file.path(OUT, "S12_kpc_ndm_pairs.csv"))
cat("\nS12: pairs", nrow(pairs), "\n")

summ <- fread(file.path(TXT, "kpc_ndm_platform_concordance_summary.csv"))
fwrite(summ, file.path(OUT, "S12_summary.csv"))
print(summ)
