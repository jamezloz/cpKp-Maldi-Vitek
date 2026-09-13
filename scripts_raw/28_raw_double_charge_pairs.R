suppressPackageStartupMessages(library(data.table))
source("scripts_raw/00_shared_utils_raw.R")

REL_TOL   <- 5e-4
MIN_R     <- 0.80
PROTON    <- 1.007

REPORTED_MZ <- c(2416, 2662, 2692, 3180, 3348, 3761, 3853, 4154, 5249, 5555,
                 5564, 5804, 7705, 11109, 11110, 11126, 11185)

ds   <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_peak_strain_treatment_dataset.rds"))
feat <- fread(file.path(RAW_MODULAR_TAB, "tophat_biomarkers", "tophat_strain_kpc_ndm.csv"))$mz
mz_all <- as.numeric(colnames(ds$matrix))
keep   <- vapply(feat, function(v) which.min(abs(mz_all - v)), integer(1))
X  <- ds$matrix[, keep]
mz <- mz_all[keep]
C  <- cor(X)
P  <- X > 0
cat(sprintf("Obs: %d | features: %d\n", nrow(X), length(mz)))

pairs <- list()
for (i in seq_along(mz)) {
  target <- 2 * mz[i] - PROTON
  js <- which(abs(mz - target) <= REL_TOL * target)
  for (j in js) {
    if (C[i, j] < MIN_R) next
    pairs[[length(pairs) + 1]] <- data.table(
      mz_z2 = round(mz[i], 1), mz_z1 = round(mz[j], 1),
      mass_approx = round(mz[j] - 1),
      pearson = round(C[i, j], 3),
      codetection = round(mean(P[, i] & P[, j]), 3))
  }
}
res <- rbindlist(pairs)

res[, dev := abs(mz_z1 - (2 * mz_z2 - PROTON))]
setorder(res, dev)
res <- res[!duplicated(mz_z2) & !duplicated(mz_z1)]
res[, dev := NULL]
res[, z2_published := vapply(mz_z2, function(v) any(abs(REPORTED_MZ - v) <= 1.5), logical(1))]
res[, z1_published := vapply(mz_z1, function(v) any(abs(REPORTED_MZ - v) <= 1.5), logical(1))]
setorder(res, -pearson)
cat(sprintf("Pares con r >= %.2f: %d (%d con algun miembro reportado)\n",
            MIN_R, nrow(res), sum(res$z2_published | res$z1_published)))
print(res[z2_published | z1_published])
save_raw_table(res, "raw_double_charge_pairs")
fwrite(res, file.path(dirname(RAW_MODULAR_OUT), "PAPER_SUPP_TABLES", "S13_double_charge_pairs.csv"))
