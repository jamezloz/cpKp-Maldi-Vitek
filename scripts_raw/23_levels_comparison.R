suppressPackageStartupMessages({
  library(data.table)
})

BASE     <- "resultados/RAW_modular"
OBJ_DIR  <- file.path(BASE, "objects")
OUT_DIR  <- file.path(BASE, "tables/levels_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREV_MIN <- 0.05

TXT_RS <- c(2416, 2662, 2692, 3180, 4154)
TXT_KN <- c(3762, 5250, 5556, 11110, 11126)

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  sum(sign(outer(x, y, "-"))) / (length(x) * length(y))
}

mat_rep <- readRDS(file.path(OBJ_DIR, "raw_peak_matrix.rds"))
det_rep <- readRDS(file.path(OBJ_DIR, "raw_peak_detected.rds"))
meta    <- as.data.table(readRDS(file.path(OBJ_DIR, "raw_metadata.rds")))
meta[, sample := basename(sample)]
stopifnot(nrow(mat_rep) == nrow(meta), identical(dim(mat_rep), dim(det_rep)))
cat(sprintf("Replicados: %d espectros, %d cepas, %d features\n",
            nrow(mat_rep), uniqueN(meta$strain_id), ncol(mat_rep)))

ds_canon  <- readRDS(file.path(OBJ_DIR, "raw_peak_strain_treatment_dataset.rds"))
res_c     <- ds_canon$metadata$label_res_sen == "RES"
prev_c    <- pmax(colMeans(ds_canon$detected[res_c, , drop = FALSE]),
                  colMeans(ds_canon$detected[!res_c, , drop = FALSE]))
FEATS <- colnames(ds_canon$matrix)[prev_c >= PREV_MIN]
cat(sprintf("Conjunto canonico de features: %d\n\n", length(FEATS)))

aggregate_level <- function(rows, keys) {
  idx_by <- split(rows, keys[rows])
  X <- do.call(rbind, lapply(idx_by, function(ix) colMeans(mat_rep[ix, , drop = FALSE])))
  D <- do.call(rbind, lapply(idx_by, function(ix) {
    as.integer(colMeans(det_rep[ix, , drop = FALSE]) >= 0.5)
  }))
  colnames(D) <- colnames(mat_rep)
  M <- rbindlist(lapply(idx_by, function(ix) meta[ix[1], .(strain_id, treatment,
                                                           mechanism, label_res_sen)]))
  list(X = X[, FEATS, drop = FALSE], D = D[, FEATS, drop = FALSE], meta = M)
}

build_level <- function(name) {
  all_rows <- seq_len(nrow(meta))
  switch(name,
    replicate = list(X = mat_rep[, FEATS, drop = FALSE],
                     D = det_rep[, FEATS, drop = FALSE] * 1L,
                     meta = meta[, .(strain_id, treatment, mechanism, label_res_sen)]),
    strain_treatment = aggregate_level(all_rows,
                         paste(meta$strain_id, meta$treatment, sep = "__")),
    strain_full      = aggregate_level(all_rows, meta$strain_id),
    col_only         = aggregate_level(which(meta$treatment == "COL"), meta$strain_id),
    ext_only         = aggregate_level(which(meta$treatment == "EXT"), meta$strain_id)
  )
}

run_task <- function(lv, task) {
  m <- lv$meta
  if (task == "res_sen") {
    keep <- m$label_res_sen %in% c("RES", "SEN"); y <- m$label_res_sen[keep]; pos <- "RES"
  } else {
    keep <- m$mechanism %in% c("KPC", "NDM");    y <- m$mechanism[keep];     pos <- "KPC"
  }
  X <- lv$X[keep, , drop = FALSE]; D <- lv$D[keep, , drop = FALSE]
  strains <- m$strain_id[keep]

  out <- rbindlist(lapply(colnames(X), function(f) {
    v <- X[, f]; b <- D[, f]
    pw <- tryCatch(wilcox.test(v[y == pos], v[y != pos], exact = FALSE)$p.value,
                   error = function(e) NA_real_)
    tab <- table(factor(b, 0:1), factor(y == pos, c(FALSE, TRUE)))
    pf <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    data.table(mz = as.numeric(f), p_wilcox = pw, p_fisher = pf,
               cliff_d = cliffs_delta(v[y == pos], v[y != pos]),
               prev_pos = mean(b[y == pos]), prev_neg = mean(b[y != pos]))
  }))
  out[, fdr_wilcox := p.adjust(p_wilcox, "fdr")]
  out[, fdr_fisher := p.adjust(p_fisher, "fdr")]
  out[, `:=`(n_obs = sum(keep), n_cepas = uniqueN(strains), task = task)]
  out[order(fdr_wilcox)]
}

LEVELS <- c("replicate", "strain_treatment", "strain_full", "col_only", "ext_only")
summ <- list(); detail <- list()

for (nm in LEVELS) {
  lv <- build_level(nm)

  ri <- lv$meta$label_res_sen == "RES"
  pv <- pmax(colMeans(lv$D[ri, , drop = FALSE]), colMeans(lv$D[!ri, , drop = FALSE]))
  n_own <- sum(pv >= PREV_MIN)

  for (task in c("res_sen", "kpc_ndm")) {
    r <- run_task(lv, task)
    r[, nivel := nm]
    detail[[paste(nm, task)]] <- r
    summ[[paste(nm, task)]] <- data.table(
      nivel = nm, task = task,
      n_obs = r$n_obs[1], n_cepas = r$n_cepas[1],
      obs_por_cepa = round(r$n_obs[1] / r$n_cepas[1], 2),
      feats_filtro_propio = n_own,
      sig_wilcox = sum(r$fdr_wilcox < 0.05, na.rm = TRUE),
      sig_fisher = sum(r$fdr_fisher < 0.05, na.rm = TRUE),
      cliff_grande = sum(abs(r$cliff_d) > 0.474, na.rm = TRUE)
    )
    cat(sprintf("%-17s %-8s n=%4d (%3d cepas, %.2f/cepa)  Wilcox=%3d  Fisher=%3d  |d|>0.474=%3d\n",
                nm, task, r$n_obs[1], r$n_cepas[1], r$n_obs[1] / r$n_cepas[1],
                sum(r$fdr_wilcox < 0.05, na.rm = TRUE),
                sum(r$fdr_fisher < 0.05, na.rm = TRUE),
                sum(abs(r$cliff_d) > 0.474, na.rm = TRUE)))
  }
}

summary_dt <- rbindlist(summ)
detail_dt  <- rbindlist(detail)
fwrite(summary_dt, file.path(OUT_DIR, "levels_summary.csv"))
fwrite(detail_dt,  file.path(OUT_DIR, "levels_detail_all_features.csv"))

cat("\n== Marcadores de referencia: FDR Wilcoxon y Cliff por nivel ==\n")
ref <- rbindlist(lapply(c("res_sen", "kpc_ndm"), function(tsk) {
  tg <- if (tsk == "res_sen") TXT_RS else TXT_KN
  rbindlist(lapply(LEVELS, function(nm) {
    d <- detail_dt[nivel == nm & task == tsk]
    rbindlist(lapply(tg, function(t) {
      near <- d[abs(mz - t) <= 12]
      if (!nrow(near)) return(NULL)
      b <- near[which.min(abs(near$mz - t))]
      data.table(task = tsk, marcador = t, nivel = nm, mz = round(b$mz, 1),
                 fdr_wilcox = signif(b$fdr_wilcox, 3), cliff_d = round(b$cliff_d, 3),
                 sig = !is.na(b$fdr_wilcox) & b$fdr_wilcox < 0.05)
    }))
  }))
}))
fwrite(ref, file.path(OUT_DIR, "levels_reference_markers.csv"))
print(dcast(ref, task + marcador ~ nivel, value.var = "sig"))

cat("\n== Concordancia de ranking (top-50 Wilcoxon) vs cepa x protocolo ==\n")
for (tsk in c("res_sen", "kpc_ndm")) {
  base_top <- head(detail_dt[nivel == "strain_treatment" & task == tsk][order(fdr_wilcox)]$mz, 50)
  for (nm in setdiff(LEVELS, "strain_treatment")) {
    o <- head(detail_dt[nivel == nm & task == tsk][order(fdr_wilcox)]$mz, 50)
    cat(sprintf("  %-8s %-17s solapamiento top-50: %d/50\n", tsk, nm,
                length(intersect(base_top, o))))
  }
}

cat("\nGuardado en:", OUT_DIR, "\n")
