suppressPackageStartupMessages({
  library(MALDIquant)
  library(MALDIquantForeign)
  library(data.table)
})

OBJ_DIR  <- "resultados/RAW_modular/objects"
OUT_DIR  <- "resultados/RAW_modular/tables/tophat_biomarkers"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREV_MIN  <- 0.05
N_BOOT    <- 200
BOOT_FRAC <- 0.50

TXT_STABLE_RS  <- c(2416, 2662, 2692, 3180, 4154)
TXT_STABLE_KN  <- c(3762, 5250, 5554, 5556, 11108, 11110, 11126)

ds <- readRDS(file.path(OBJ_DIR, "raw_peak_strain_treatment_dataset.rds"))
meta_strain <- as.data.table(ds$metadata)
X_all   <- ds$matrix
det_all <- ds$detected
stopifnot(nrow(X_all) == nrow(meta_strain),
          identical(dim(X_all), dim(det_all)),
          all(det_all %in% c(0, 1)))
cat(sprintf("Dataset canonico: %d obs (cepa x protocolo), %d cepas, %d features\n",
            nrow(X_all), uniqueN(meta_strain$strain_id), ncol(X_all)))

res_i     <- meta_strain$label_res_sen == "RES"
prev_max  <- pmax(colMeans(det_all[res_i, , drop = FALSE]),
                  colMeans(det_all[!res_i, , drop = FALSE]))
keep_prev <- prev_max >= PREV_MIN

mat_strain     <- X_all[,   keep_prev, drop = FALSE]
mat_bin_strain <- det_all[, keep_prev, drop = FALSE]
cat(sprintf("  Prevalencia maxima por clase >= %.0f%%: %d de %d features\n",
            PREV_MIN * 100, sum(keep_prev), length(prev_max)))

cat("  KPC:", sum(meta_strain$mechanism=="KPC"),
    "| NDM:", sum(meta_strain$mechanism=="NDM"),
    "| RES:", sum(meta_strain$label_res_sen=="RES"),
    "| SEN:", sum(meta_strain$label_res_sen=="SEN"), "\n")

cat("  Binary prevalence check (majority-vote, ≥50% reps):\n")
for (tgt in c(5249, 5556, 11109, 3762)) {
  near <- which(abs(as.numeric(colnames(mat_bin_strain)) - tgt) <= 10)
  if (length(near) > 0) {
    idx_c <- near[which.min(abs(as.numeric(colnames(mat_bin_strain))[near] - tgt))]
    kpc_p <- mean(mat_bin_strain[meta_strain$mechanism=="KPC", idx_c])
    ndm_p <- mean(mat_bin_strain[meta_strain$mechanism=="NDM", idx_c])
    cat(sprintf("    ~%d Da: KPC=%.1f%% NDM=%.1f%%\n", tgt, kpc_p*100, ndm_p*100))
  }
}

cliffs_delta <- function(x, y) {

  nx <- length(x); ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  mat <- outer(x, y, function(a, b) sign(a - b))
  sum(mat) / (nx * ny)
}

run_wilcoxon_strain <- function(mat, y, label_pos, log2_transform = FALSE) {
  if (log2_transform) mat <- log2(mat + 1e-10)
  rbindlist(lapply(colnames(mat), function(cn) {
    v <- mat[, cn]
    w <- tryCatch(wilcox.test(v[y == label_pos], v[y != label_pos], exact = FALSE),
                  error = function(e) NULL)
    if (is.null(w)) return(NULL)
    cd <- cliffs_delta(v[y == label_pos], v[y != label_pos])
    data.table(
      mz         = as.numeric(cn),
      p_wilcox   = w$p.value,
      median_pos = median(v[y == label_pos]),
      median_neg = median(v[y != label_pos]),
      cliff_d    = cd,
      n_pos      = sum(y == label_pos),
      n_neg      = sum(y != label_pos)
    )
  }))
}

run_fisher_strain <- function(mat_bin, y, label_pos) {
  rbindlist(lapply(colnames(mat_bin), function(cn) {
    b   <- mat_bin[, cn]

    tab <- table(factor(b, 0:1), factor(y == label_pos, c(FALSE, TRUE)))
    ft  <- tryCatch(fisher.test(tab, simulate.p.value = FALSE), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    data.table(
      mz      = as.numeric(cn),
      p_fisher = ft$p.value,
      OR       = as.numeric(ft$estimate),
      pct_pos  = mean(b[y == label_pos]),
      pct_neg  = mean(b[y != label_pos])
    )
  }))
}

combine_tests <- function(wlx, fsh, label_pos, label_neg) {
  m <- merge(wlx, fsh, by = "mz", all = TRUE)
  m[, fdr_wilcox := p.adjust(p_wilcox, "fdr")]
  m[, fdr_fisher := p.adjust(p_fisher,  "fdr")]
  m[, sig_wilcox := !is.na(fdr_wilcox) & fdr_wilcox < 0.05]
  m[, sig_fisher := !is.na(fdr_fisher) & fdr_fisher  < 0.05]
  m[, sig_both   := sig_wilcox & sig_fisher]
  m[, log2fc     := log2((median_pos + 1e-10) / (median_neg + 1e-10))]
  m[, combined_score := ifelse(is.na(p_wilcox), 0, -log10(p_wilcox)) +
                         ifelse(is.na(p_fisher),  0, -log10(p_fisher))]
  m[, group_label := paste0(label_pos, " vs ", label_neg)]
  m[order(-combined_score)]
}

check_txt_overlap <- function(mz_vec, targets, tol = 10) {
  sapply(targets, function(t) any(abs(mz_vec - t) <= tol))
}

bootstrap_stability <- function(mat, y, strain, label_pos,
                                n_boot = N_BOOT, frac = BOOT_FRAC) {
  feat_names <- colnames(mat)
  sel_count  <- integer(length(feat_names))
  names(sel_count) <- feat_names

  strains_pos <- unique(strain[y == label_pos])
  strains_neg <- unique(strain[y != label_pos])

  for (b in seq_len(n_boot)) {
    set.seed(b)
    s_pos  <- sample(strains_pos, max(2L, ceiling(length(strains_pos) * frac)))
    s_neg  <- sample(strains_neg, max(2L, ceiling(length(strains_neg) * frac)))
    si     <- which(strain %in% c(s_pos, s_neg))
    yi     <- y[si]
    mi     <- mat[si, , drop = FALSE]

    pvals <- vapply(colnames(mi), function(cn) {
      v <- mi[, cn]
      tryCatch(wilcox.test(v[yi==label_pos], v[yi!=label_pos], exact=FALSE)$p.value,
               error = function(e) 1.0)
    }, numeric(1))

    fdr <- p.adjust(pvals, "fdr")
    sig <- names(which(!is.na(fdr) & fdr < 0.05))
    sel_count[sig] <- sel_count[sig] + 1L
  }

  data.table(
    mz = as.numeric(feat_names),
    feature = feat_names,
    stability = sel_count / n_boot
  )[order(-stability)]
}

cat("\n══════════════════════════════════════\n")
cat("KPC vs NDM (cepa x protocolo)\n")
cat("══════════════════════════════════════\n")

idx_kn    <- meta_strain$mechanism %in% c("KPC", "NDM")
y_kn      <- meta_strain$mechanism[idx_kn]
m_kn      <- mat_strain[idx_kn, , drop = FALSE]
mb_kn     <- mat_bin_strain[idx_kn, , drop = FALSE]
strain_kn <- meta_strain$strain_id[idx_kn]

wlx_kn   <- run_wilcoxon_strain(m_kn, y_kn, "KPC", log2_transform = FALSE)
wlx_kn_log <- run_wilcoxon_strain(m_kn, y_kn, "KPC", log2_transform = TRUE)
fsh_kn   <- run_fisher_strain(mb_kn, y_kn, "KPC")

wlx_kn[, p_wilcox_log2 := wlx_kn_log$p_wilcox[match(mz, wlx_kn_log$mz)]]
wlx_kn[, cliff_d_log2  := wlx_kn_log$cliff_d[match(mz, wlx_kn_log$mz)]]

bio_kn <- combine_tests(wlx_kn, fsh_kn, "KPC", "NDM")

cat(sprintf("  Wilcoxon sig (FDR<0.05): %d\n", sum(bio_kn$sig_wilcox, na.rm=TRUE)))
cat(sprintf("  Fisher sig  (FDR<0.05): %d\n", sum(bio_kn$sig_fisher,  na.rm=TRUE)))
cat(sprintf("  Both sig:               %d\n", sum(bio_kn$sig_both,    na.rm=TRUE)))
cat(sprintf("  |Cliff d|>0.33 (medium): %d\n", sum(abs(bio_kn$cliff_d)>0.33, na.rm=TRUE)))
cat(sprintf("  |Cliff d|>0.474 (large): %d\n", sum(abs(bio_kn$cliff_d)>0.474, na.rm=TRUE)))

ovl_kn_wlx  <- check_txt_overlap(bio_kn[sig_wilcox==TRUE]$mz, TXT_STABLE_KN)
ovl_kn_both <- check_txt_overlap(bio_kn[sig_both==TRUE]$mz,    TXT_STABLE_KN)
cat(sprintf("  TXT overlap (Wilcoxon): %d/%d | (Both): %d/%d\n",
            sum(ovl_kn_wlx), length(TXT_STABLE_KN),
            sum(ovl_kn_both), length(TXT_STABLE_KN)))

cat("\n  Top 15 by combined score (sig_both):\n")
top_kn <- bio_kn[sig_both == TRUE][order(-combined_score)][1:min(15, .N)]
cat(sprintf("  %-8s  KPC%%  NDM%%   OR    Cliff_d  log2FC  p_wilcox   p_fisher\n", "mz"))
for (i in seq_len(nrow(top_kn))) {
  cat(sprintf("  %-8.1f  %4.1f%%  %4.1f%%  %5.1f  %+6.3f  %+6.2f  %.2e   %.2e\n",
    top_kn$mz[i],
    ifelse(is.na(top_kn$pct_pos[i]),0,top_kn$pct_pos[i]*100),
    ifelse(is.na(top_kn$pct_neg[i]),0,top_kn$pct_neg[i]*100),
    ifelse(is.na(top_kn$OR[i]),0,top_kn$OR[i]),
    ifelse(is.na(top_kn$cliff_d[i]),0,top_kn$cliff_d[i]),
    ifelse(is.na(top_kn$log2fc[i]),0,top_kn$log2fc[i]),
    top_kn$p_wilcox[i], top_kn$p_fisher[i]))
}

cat("\n  Running bootstrap stability (200 iter, bloqueado por cepa)...\n")
stab_kn <- bootstrap_stability(m_kn, y_kn, strain_kn, "KPC")
stab_kn[, txt_target := sapply(mz, function(m) any(abs(TXT_STABLE_KN - m) <= 10))]
cat(sprintf("  Stability=1.00: %d features\n", sum(stab_kn$stability == 1.0)))
cat("  TXT targets stability:\n")
for (t in TXT_STABLE_KN) {
  near <- stab_kn[abs(mz - t) <= 10][1]
  cat(sprintf("    %5d Da → found: %6.1f Da, stab=%.3f\n",
              t, ifelse(is.na(near$mz[1]),t,near$mz[1]),
              ifelse(is.na(near$stability[1]),0,near$stability[1])))
}

bio_kn <- merge(bio_kn, stab_kn[, .(mz=as.numeric(feature), stability)], by="mz", all.x=TRUE)
fwrite(bio_kn[order(-combined_score)], file.path(OUT_DIR, "tophat_strain_kpc_ndm.csv"))
fwrite(stab_kn, file.path(OUT_DIR, "bootstrap_stability_tophat_kpc_ndm.csv"))

cat("\n══════════════════════════════════════\n")
cat("RES vs SEN (cepa x protocolo)\n")
cat("══════════════════════════════════════\n")

idx_rs    <- meta_strain$label_res_sen %in% c("RES", "SEN")
y_rs      <- meta_strain$label_res_sen[idx_rs]
m_rs      <- mat_strain[idx_rs, , drop = FALSE]
mb_rs     <- mat_bin_strain[idx_rs, , drop = FALSE]
strain_rs <- meta_strain$strain_id[idx_rs]

wlx_rs     <- run_wilcoxon_strain(m_rs, y_rs, "RES", log2_transform = FALSE)
wlx_rs_log <- run_wilcoxon_strain(m_rs, y_rs, "RES", log2_transform = TRUE)
fsh_rs     <- run_fisher_strain(mb_rs, y_rs, "RES")

wlx_rs[, p_wilcox_log2 := wlx_rs_log$p_wilcox[match(mz, wlx_rs_log$mz)]]
wlx_rs[, cliff_d_log2  := wlx_rs_log$cliff_d[match(mz, wlx_rs_log$mz)]]

bio_rs <- combine_tests(wlx_rs, fsh_rs, "RES", "SEN")

cat(sprintf("  Wilcoxon sig (FDR<0.05): %d\n", sum(bio_rs$sig_wilcox, na.rm=TRUE)))
cat(sprintf("  Fisher sig  (FDR<0.05): %d\n", sum(bio_rs$sig_fisher,  na.rm=TRUE)))
cat(sprintf("  Both sig:               %d\n", sum(bio_rs$sig_both,    na.rm=TRUE)))
cat(sprintf("  |Cliff d|>0.33 (medium): %d\n", sum(abs(bio_rs$cliff_d)>0.33, na.rm=TRUE)))
cat(sprintf("  |Cliff d|>0.474 (large): %d\n", sum(abs(bio_rs$cliff_d)>0.474, na.rm=TRUE)))

ovl_rs_wlx  <- check_txt_overlap(bio_rs[sig_wilcox==TRUE]$mz, TXT_STABLE_RS)
ovl_rs_both <- check_txt_overlap(bio_rs[sig_both==TRUE]$mz,    TXT_STABLE_RS)
cat(sprintf("  TXT overlap (Wilcoxon): %d/%d | (Both): %d/%d\n",
            sum(ovl_rs_wlx), length(TXT_STABLE_RS),
            sum(ovl_rs_both), length(TXT_STABLE_RS)))

cat("\n  Top 15 by combined score (sig_both):\n")
top_rs <- bio_rs[sig_both == TRUE][order(-combined_score)][1:min(15, .N)]
cat(sprintf("  %-8s  RES%%  SEN%%   OR    Cliff_d  log2FC  p_wilcox   p_fisher\n", "mz"))
for (i in seq_len(nrow(top_rs))) {
  cat(sprintf("  %-8.1f  %4.1f%%  %4.1f%%  %5.1f  %+6.3f  %+6.2f  %.2e   %.2e\n",
    top_rs$mz[i],
    ifelse(is.na(top_rs$pct_pos[i]),0,top_rs$pct_pos[i]*100),
    ifelse(is.na(top_rs$pct_neg[i]),0,top_rs$pct_neg[i]*100),
    ifelse(is.na(top_rs$OR[i]),0,top_rs$OR[i]),
    ifelse(is.na(top_rs$cliff_d[i]),0,top_rs$cliff_d[i]),
    ifelse(is.na(top_rs$log2fc[i]),0,top_rs$log2fc[i]),
    top_rs$p_wilcox[i], top_rs$p_fisher[i]))
}

cat("\n  Running bootstrap stability (200 iter, bloqueado por cepa)...\n")
stab_rs <- bootstrap_stability(m_rs, y_rs, strain_rs, "RES")
stab_rs[, txt_target := sapply(mz, function(m) any(abs(TXT_STABLE_RS - m) <= 10))]
cat(sprintf("  Stability=1.00: %d features\n", sum(stab_rs$stability == 1.0)))
cat("  TXT targets stability:\n")
for (t in TXT_STABLE_RS) {
  near <- stab_rs[abs(mz - t) <= 10][1]
  cat(sprintf("    %5d Da → found: %6.1f Da, stab=%.3f\n",
              t, ifelse(is.na(near$mz[1]),t,near$mz[1]),
              ifelse(is.na(near$stability[1]),0,near$stability[1])))
}

bio_rs <- merge(bio_rs, stab_rs[, .(mz=as.numeric(feature), stability)], by="mz", all.x=TRUE)
fwrite(bio_rs[order(-combined_score)], file.path(OUT_DIR, "tophat_strain_res_sen.csv"))
fwrite(stab_rs, file.path(OUT_DIR, "bootstrap_stability_tophat_res_sen.csv"))

cat("\n══════════════════════════════════════════════════════════\n")
cat("PAPER TABLE: Top biomarkers sig_both + |Cliff|>0.33\n")
cat("══════════════════════════════════════════════════════════\n")

make_paper_table <- function(bio, stab, group_label, txt_targets) {
  dt <- bio[sig_both == TRUE & abs(cliff_d) > 0.33][order(-combined_score)]
  if (nrow(dt) == 0) return(NULL)
  dt[, stability := stab$stability[match(as.character(round(mz,1)),
                                          as.character(round(as.numeric(stab$feature),1)))]]
  dt[, txt_confirmed := sapply(mz, function(m) any(abs(txt_targets - m) <= 10))]
  dt[, group := group_label]
  dt[, .(group, mz=round(mz,1), pct_pos=round(pct_pos*100,1), pct_neg=round(pct_neg*100,1),
         OR=round(OR,2), cliff_d=round(cliff_d,3), log2fc=round(log2fc,2),
         fdr_wilcox=signif(fdr_wilcox,3), fdr_fisher=signif(fdr_fisher,3),
         stability=round(stability,2), txt_confirmed)]
}

tbl_kn <- make_paper_table(bio_kn, stab_kn, "KPC_vs_NDM", TXT_STABLE_KN)
tbl_rs <- make_paper_table(bio_rs, stab_rs, "RES_vs_SEN", TXT_STABLE_RS)
tbl_paper <- rbindlist(list(tbl_kn, tbl_rs), fill = TRUE)

print(tbl_paper)
fwrite(tbl_paper, file.path(OUT_DIR, "paper_biomarker_table_tophat_strain.csv"))

cat("\n\n══════════════════════════════════════════════════════════\n")

cat("RESUMEN (pipeline unico, cepa x protocolo)\n")
cat("══════════════════════════════════════════════════════════\n")
cat("                 KPC/NDM                  RES/SEN\n")
cat("                 Wilcox  Fisher  Both     Wilcox  Fisher  Both\n")
cat(sprintf("Significativos:  %-6d  %-6d  %-6d   %-6d  %-6d  %-6d\n",
    sum(bio_kn$sig_wilcox,na.rm=T), sum(bio_kn$sig_fisher,na.rm=T), sum(bio_kn$sig_both,na.rm=T),
    sum(bio_rs$sig_wilcox,na.rm=T), sum(bio_rs$sig_fisher,na.rm=T), sum(bio_rs$sig_both,na.rm=T)))

cat("\nAll results saved to:", OUT_DIR, "\n")
