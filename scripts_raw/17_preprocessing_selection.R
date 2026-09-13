source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(data.table)
  library(MALDIquant)
})

SMOOTH_CACHE <- file.path(RAW_MODULAR_OBJ, "raw_specs_smoothed_only.rds")
OUT_DIR <- file.path(RAW_MODULAR_TAB, "preprocessing_selection")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_FREQ <- 0.5

baselines <- list(
  snip25  = function(s) removeBaseline(s, method = "SNIP", iterations = 25),
  snip50  = function(s) removeBaseline(s, method = "SNIP", iterations = 50),
  snip100 = function(s) removeBaseline(s, method = "SNIP", iterations = 100),
  tophat  = function(s) removeBaseline(s, method = "TopHat")
)

repro_stats <- function(detected, groups) {
  idx <- split(seq_len(nrow(detected)), groups)
  idx <- idx[vapply(idx, length, integer(1)) >= 2]
  if (!length(idx)) return(NULL)
  per_group <- lapply(idx, function(i) {
    sub <- detected[i, , drop = FALSE]
    presente <- colSums(sub) > 0
    if (!any(presente)) return(NULL)
    frac <- colMeans(sub[, presente, drop = FALSE])
    data.table(n_detectados = sum(presente),
               n_consistentes = sum(frac >= MIN_FREQ),
               consistencia = mean(frac))
  })
  rbindlist(per_group[!vapply(per_group, is.null, logical(1))])
}

main <- function() {
  set.seed(42)
  stopifnot(file.exists(SMOOTH_CACHE))
  specs_smooth <- readRDS(SMOOTH_CACHE)
  meta <- as.data.table(readRDS(file.path(RAW_MODULAR_OBJ, "raw_metadata.rds")))
  stopifnot(length(specs_smooth) == nrow(meta),
            identical(names(specs_smooth), meta$sample))

  groups <- paste(meta$strain_id, meta$treatment, sep = "__")
  n_con_replicas <- sum(table(groups) >= 2)
  cat("espectros:", nrow(meta), "| grupos cepa x protocolo con >=2 replicas:",
      n_con_replicas, "\n\n")

  out <- list()
  for (nm in names(baselines)) {
    cat("===", nm, "===\n")
    specs <- lapply(specs_smooth, baselines[[nm]])
    peaks <- lapply(specs, function(s) detectPeaks(s, SNR = estimate_snr(s), halfWindowSize = 10))
    snr_med <- median(unlist(lapply(peaks, function(p) MALDIquant::snr(p))), na.rm = TRUE)

    pb <- binPeaks(peaks, tolerance = 0.003)

    det <- !is.na(intensityMatrix(pb))
    rownames(det) <- names(specs)
    det <- det[meta$sample, , drop = FALSE]

    r <- repro_stats(det, groups)
    if (is.null(r)) { cat("  sin grupos evaluables\n"); next }
    out[[length(out) + 1]] <- data.table(
      preprocesado = nm,
      n_features_total = ncol(det),
      picos_por_espectro = round(mean(rowSums(det)), 1),
      snr_mediano = round(snr_med, 2),
      consistencia = round(mean(r$consistencia), 4),
      n_consistentes_medio = round(mean(r$n_consistentes), 1),
      pct_consistentes = round(mean(r$n_consistentes / r$n_detectados), 4))
    cat("  features:", ncol(det), "| consistencia:", round(mean(r$consistencia), 3),
        "| picos consistentes:", round(mean(r$n_consistentes), 1),
        "| SNR:", round(snr_med, 1), "\n")
  }

  res <- rbindlist(out)
  setorder(res, -consistencia)
  fwrite(res, file.path(OUT_DIR, "preprocesado_reproducibilidad.csv"))
  cat("\n=== RESULTADO (ordenado por reproducibilidad) ===\n")
  print(res)
  cat("\nconsistencia   = fraccion media de replicas en que aparece cada pico\n")
  cat("n_consistentes = picos presentes en >=", round(100*MIN_FREQ), "% de las replicas\n")
  cat("\nLas tres metricas se leen juntas: mucha consistencia con pocos picos no\n")
  cat("sirve, y muchos picos poco consistentes tampoco.\n")
  cat("\nSalidas en:", OUT_DIR, "\n")
}

main()
