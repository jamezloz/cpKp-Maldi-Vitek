suppressPackageStartupMessages(library(data.table))

BASE    <- "resultados/RAW_modular"
OBJ_DIR <- file.path(BASE, "objects")
OUT_DIR <- file.path(BASE, "tables/panel_complementarity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_SPEC  <- 0.98
MAX_PANEL <- 6
N_BOOT    <- 200

ds   <- readRDS(file.path(OBJ_DIR, "raw_peak_strain_treatment_dataset.rds"))
det  <- ds$detected
meta <- as.data.table(ds$metadata)
mz   <- as.numeric(colnames(det))

keep <- meta$mechanism %in% c("KPC", "NDM")
det  <- det[keep, , drop = FALSE]
meta <- meta[keep]
cat(sprintf("Obs cepa x protocolo: %d (%d cepas) | KPC=%d NDM=%d\n",
            nrow(det), uniqueN(meta$strain_id),
            sum(meta$mechanism == "KPC"), sum(meta$mechanism == "NDM")))

prev_kpc <- colMeans(det[meta$mechanism == "KPC", , drop = FALSE])
prev_ndm <- colMeans(det[meta$mechanism == "NDM", , drop = FALSE])
cand <- which(pmax(prev_kpc, prev_ndm) >= 0.05)
cat(sprintf("Features candidatas (prevalencia maxima por clase >= 5%%): %d\n\n", length(cand)))

greedy_panel <- function(D, y, target, other, cand_idx, min_spec = MIN_SPEC,
                         max_k = MAX_PANEL, seed_feat = NULL) {
  is_t <- y == target; is_o <- y == other
  panel <- integer(0)
  hit   <- rep(FALSE, nrow(D))
  steps <- list()

  if (!is.null(seed_feat)) {
    panel <- seed_feat
    hit <- D[, seed_feat] == 1
    steps[[1]] <- data.table(paso = 0L, mz = round(mz[seed_feat], 1),
                             sens = mean(hit[is_t]), spec = 1 - mean(hit[is_o]),
                             nuevos = sum(hit[is_t]), semilla = TRUE)
  }

  for (k in seq_len(max_k - length(panel))) {
    best <- NULL
    for (j in setdiff(cand_idx, panel)) {
      h2 <- hit | (D[, j] == 1)
      sp <- 1 - mean(h2[is_o])
      if (sp < min_spec) next
      se <- mean(h2[is_t])
      nuevos <- sum(h2[is_t]) - sum(hit[is_t])
      if (nuevos <= 0) next
      if (is.null(best) || se > best$se) best <- list(j = j, se = se, sp = sp, nuevos = nuevos)
    }
    if (is.null(best)) break
    panel <- c(panel, best$j)
    hit   <- hit | (D[, best$j] == 1)
    steps[[length(steps) + 1]] <- data.table(
      paso = length(panel) - as.integer(!is.null(seed_feat)),
      mz = round(mz[best$j], 1), sens = best$se, spec = best$sp,
      nuevos = best$nuevos, semilla = FALSE)
  }
  list(panel = panel, steps = rbindlist(steps, fill = TRUE), hit = hit)
}

eval_panel <- function(D, y, panel, target, other) {
  h <- if (length(panel) == 1) D[, panel] == 1 else rowSums(D[, panel, drop = FALSE]) > 0
  c(sens = mean(h[y == target]), spec = 1 - mean(h[y == other]))
}

i_11109 <- which.min(abs(mz - 11109))
i_5555  <- which.min(abs(mz - 5555))
cat("=========== KPC: complementariedad sobre 11109 Da ===========\n")
cat(sprintf("Semilla: %.1f Da  sens=%.3f  spec=%.3f\n",
            mz[i_11109], mean(det[meta$mechanism=="KPC", i_11109]),
            1 - mean(det[meta$mechanism=="NDM", i_11109])))
u <- (det[, i_11109] == 1) | (det[, i_5555] == 1)
cat(sprintf("11109 O 5555 (dos cargas del mismo ente): sens=%.3f spec=%.3f\n\n",
            mean(u[meta$mechanism=="KPC"]), 1 - mean(u[meta$mechanism=="NDM"])))

g_kpc <- greedy_panel(det, meta$mechanism, "KPC", "NDM", cand, seed_feat = i_11109)
print(g_kpc$steps)
cat(sprintf("\nPanel KPC final (%d marcadores): %s\n", length(g_kpc$panel),
            paste(round(mz[g_kpc$panel], 1), collapse = ", ")))
fwrite(g_kpc$steps, file.path(OUT_DIR, "kpc_greedy_steps.csv"))

cat("\n=========== NDM: mejor marcador individual ===========\n")

single <- rbindlist(lapply(cand, function(j) {
  s <- eval_panel(det, meta$mechanism, j, "NDM", "KPC")
  data.table(mz = round(mz[j], 1), sens = s["sens"], spec = s["spec"],
             youden = s["sens"] + s["spec"] - 1)
}))
setorder(single, -youden)
print(head(single, 10))
fwrite(single, file.path(OUT_DIR, "ndm_single_markers.csv"))

cat("\n=========== NDM: panel greedy (spec >= 0.90) ===========\n")

g_ndm <- greedy_panel(det, meta$mechanism, "NDM", "KPC", cand, min_spec = 0.90)
print(g_ndm$steps)
cat(sprintf("\nPanel NDM final (%d marcadores): %s\n", length(g_ndm$panel),
            paste(round(mz[g_ndm$panel], 1), collapse = ", ")))
fwrite(g_ndm$steps, file.path(OUT_DIR, "ndm_greedy_steps.csv"))

cat("\n=========== Validacion bootstrap (seleccion dentro, eval OOB) ===========\n")
strains <- unique(meta$strain_id)
boot_eval <- function(target, other, min_spec, seed_feat = NULL) {
  res <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
    set.seed(b)
    s_in  <- sample(strains, length(strains), replace = TRUE)
    ri    <- which(meta$strain_id %in% s_in)
    oob   <- which(!meta$strain_id %in% s_in)
    if (!length(oob)) return(NULL)
    yi <- meta$mechanism[ri]; yo <- meta$mechanism[oob]
    if (uniqueN(yi) < 2 || uniqueN(yo) < 2) return(NULL)
    pk <- colMeans(det[ri[yi == "KPC"], , drop = FALSE])
    pn <- colMeans(det[ri[yi == "NDM"], , drop = FALSE])
    ci <- which(pmax(pk, pn) >= 0.05)
    if (!length(ci)) return(NULL)
    g <- greedy_panel(det[ri, , drop = FALSE], yi, target, other, ci,
                      min_spec = min_spec, seed_feat = seed_feat)
    if (!length(g$panel)) return(NULL)
    e <- eval_panel(det[oob, , drop = FALSE], yo, g$panel, target, other)
    data.table(k = length(g$panel), sens = e["sens"], spec = e["spec"])
  }))
  res
}
for (cfg in list(list(t="KPC", o="NDM", ms=MIN_SPEC, sf=i_11109, lab="KPC (semilla 11109)"),
                 list(t="NDM", o="KPC", ms=0.90,     sf=NULL,    lab="NDM (spec>=0.90)"))) {
  r <- boot_eval(cfg$t, cfg$o, cfg$ms, cfg$sf)
  if (!is.null(r) && nrow(r)) {
    cat(sprintf("%-22s  n=%d iter | sens OOB %.3f (IC95 %.3f-%.3f) | spec OOB %.3f (IC95 %.3f-%.3f) | k medio %.1f\n",
        cfg$lab, nrow(r), mean(r$sens), quantile(r$sens,.025), quantile(r$sens,.975),
        mean(r$spec), quantile(r$spec,.025), quantile(r$spec,.975), mean(r$k)))
    fwrite(r, file.path(OUT_DIR, sprintf("bootstrap_oob_%s.csv", cfg$t)))
  }
}

cat("\nGuardado en:", OUT_DIR, "\n")
