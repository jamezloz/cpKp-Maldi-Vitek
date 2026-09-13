suppressPackageStartupMessages({ library(data.table); library(MALDIquant) })

BASE    <- "resultados/RAW_modular"
OBJ_DIR <- file.path(BASE, "objects")
OUT_DIR <- file.path(BASE, "tables/improved_peaklist")
N_BOOT  <- 200
SET_PCT <- 0.5

TXT_RS <- c(2416, 2662, 2692, 3180, 4154)
TXT_KN <- c(3762, 5250, 5554, 5556, 11108, 11110, 11126)

cat("Cargando dataset canonico (cepa x protocolo)...\n")
ds   <- readRDS(file.path(OBJ_DIR, "raw_peak_strain_treatment_dataset.rds"))
meta <- as.data.table(ds$metadata)
det  <- ds$detected
stopifnot(nrow(ds$matrix) == nrow(meta), identical(dim(ds$matrix), dim(det)))
cat(sprintf("  %d obs (cepa x protocolo), %d cepas, %d features\n",
            nrow(ds$matrix), uniqueN(meta$strain_id), ncol(ds$matrix)))

res_i    <- meta$label_res_sen == "RES"
prev_max <- pmax(colMeans(det[res_i, , drop = FALSE]),
                 colMeans(det[!res_i, , drop = FALSE]))
mat <- ds$matrix[, prev_max >= 0.05, drop = FALSE]
cat("  Features tras prevalencia maxima por clase >= 5%:", ncol(mat), "\n\n")

run_boot <- function(mat, y, strain, label_pos, n_boot, set_pct) {
  levels_y <- unique(y)
  label_neg <- setdiff(levels_y, label_pos)
  sig_count <- integer(ncol(mat))
  names(sig_count) <- colnames(mat)

  strains_pos <- unique(strain[y == label_pos])
  strains_neg <- unique(strain[y == label_neg])

  for (b in seq_len(n_boot)) {
    set.seed(b)
    s_pos   <- sample(strains_pos, ceiling(length(strains_pos) * set_pct))
    s_neg   <- sample(strains_neg, ceiling(length(strains_neg) * set_pct))
    idx_s   <- which(strain %in% c(s_pos, s_neg))
    ms      <- mat[idx_s, , drop = FALSE]
    ys      <- y[idx_s]

    pvals <- sapply(seq_len(ncol(ms)), function(j) {
      v <- ms[, j]
      tryCatch(wilcox.test(v ~ factor(ys), exact = FALSE)$p.value, error = function(e) 1)
    })
    fdr   <- p.adjust(pvals, "fdr")
    sig_count <- sig_count + (!is.na(fdr) & fdr < 0.05)
  }
  data.table(mz = as.numeric(names(sig_count)),
             stability = sig_count / n_boot)
}

cat("Bootstrap RES/SEN (", N_BOOT, "iterations)...\n")
idx_rs <- meta$label_res_sen %in% c("RES","SEN")
y_rs   <- meta$label_res_sen[idx_rs]
m_rs   <- mat[idx_rs, , drop = FALSE]
strain_rs <- meta$strain_id[idx_rs]
stab_rs <- run_boot(m_rs, y_rs, strain_rs, "RES", N_BOOT, SET_PCT)
stab_rs[, in_txt := sapply(mz, function(x) any(abs(x - TXT_RS) <= 10))]
setorder(stab_rs, -stability)

cat("\nTop 20 RES/SEN stable features (stab >= 0.6):\n")
print(stab_rs[stability >= 0.6][1:min(20,.N)])
cat("TXT targets stability:\n")
for (t in TXT_RS) {
  r <- stab_rs[abs(mz-t)<=10][which.min(abs(mz-t))]
  if (nrow(r) > 0) cat(sprintf("  TXT %d -> RAW %.1f: stab=%.3f\n", t, r$mz, r$stability))
}

cat("\nBootstrap KPC/NDM (", N_BOOT, "iterations)...\n")
idx_kn <- meta$mechanism %in% c("KPC","NDM")
y_kn   <- meta$mechanism[idx_kn]
m_kn   <- mat[idx_kn, , drop = FALSE]
strain_kn <- meta$strain_id[idx_kn]
stab_kn <- run_boot(m_kn, y_kn, strain_kn, "KPC", N_BOOT, SET_PCT)
stab_kn[, in_txt := sapply(mz, function(x) any(abs(x - TXT_KN) <= 10))]
setorder(stab_kn, -stability)

cat("\nTop 20 KPC/NDM stable features (stab >= 0.6):\n")
print(stab_kn[stability >= 0.6][1:min(20,.N)])
cat("TXT targets stability:\n")
for (t in TXT_KN) {
  r <- stab_kn[abs(mz-t)<=10][which.min(abs(mz-t))]
  if (nrow(r) > 0) cat(sprintf("  TXT %d -> RAW %.1f: stab=%.3f\n", t, r$mz, r$stability))
}

fwrite(stab_rs, file.path(OUT_DIR, "bootstrap_stability_res_sen.csv"))
fwrite(stab_kn, file.path(OUT_DIR, "bootstrap_stability_kpc_ndm.csv"))
cat("\nSaved to", OUT_DIR, "\n")
