source('scripts_raw/00_shared_utils_raw.R')

suppressPackageStartupMessages({
  library(data.table)
})

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  sum(sign(outer(x, y, "-"))) / (length(x) * length(y))
}

main <- function() {
  ds <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_strain_treatment_dataset.rds'))
  X <- ds$matrix
  meta <- as.data.table(ds$metadata)
  cat(sprintf("  Nivel de analisis: cepa x protocolo (n = %d, cepas = %d)\n", nrow(X), uniqueN(meta$strain_id)))
  feat_names <- colnames(X)

  det <- ds$detected
  res_i <- meta$label_res_sen == 'RES'
  prev_max <- pmax(colMeans(det[res_i, , drop = FALSE]),
                   colMeans(det[!res_i, , drop = FALSE]))
  keep_prev <- prev_max >= 0.05
  Xg <- X[, keep_prev, drop = FALSE]
  cat(sprintf("  Prevalencia maxima por clase >= 5%%: %d de %d features\n",
              sum(keep_prev), length(prev_max)))

  keep4 <- meta$mechanism %in% c('KPC', 'NDM', 'SEN', 'KPC-NDM')
  y4 <- factor(meta$mechanism[keep4], levels = c('SEN', 'KPC', 'NDM', 'KPC-NDM'))
  kr <- lapply(colnames(Xg), function(f) data.table(feature_chr = f, p_value = kruskal.test(Xg[keep4, f] ~ y4)$p.value))
  kr <- rbindlist(kr)[order(p_value)]
  kr[, fdr := p.adjust(p_value, 'fdr')]
  save_raw_table(kr, 'raw_classic_biomarkers_kruskal_global')

  keep_res <- !is.na(meta$label_res_sen)
  y_res <- factor(meta$label_res_sen[keep_res], levels = c('SEN', 'RES'))
  res_rows <- lapply(colnames(Xg), function(f) {
    p <- wilcox.test(Xg[keep_res, f] ~ y_res)$p.value
    vals <- Xg[keep_res, f]
    med_res <- median(vals[y_res == 'RES'], na.rm = TRUE)
    med_sen <- median(vals[y_res == 'SEN'], na.rm = TRUE)
    cd <- cliffs_delta(vals[y_res == 'RES'], vals[y_res == 'SEN'])
    data.table(feature_chr = f, p_value = p, median_res = med_res, median_sen = med_sen,
               cliffs_delta = cd, direction = ifelse(med_res > med_sen, 'RES', 'SEN'))
  })
  res_dt <- rbindlist(res_rows)[order(p_value)]
  res_dt[, fdr := p.adjust(p_value, 'fdr')]
  save_raw_table(res_dt, 'raw_classic_final_biomarkers_res_sen')

  keep_mech <- meta$mechanism %in% c('KPC', 'NDM')
  y_mech <- factor(meta$mechanism[keep_mech], levels = c('NDM', 'KPC'))
  mech_rows <- lapply(colnames(Xg), function(f) {
    p <- wilcox.test(Xg[keep_mech, f] ~ y_mech)$p.value
    vals <- Xg[keep_mech, f]
    med_kpc <- median(vals[y_mech == 'KPC'], na.rm = TRUE)
    med_ndm <- median(vals[y_mech == 'NDM'], na.rm = TRUE)
    cd <- cliffs_delta(vals[y_mech == 'KPC'], vals[y_mech == 'NDM'])
    data.table(feature_chr = f, p_value = p, median_kpc = med_kpc, median_ndm = med_ndm,
               cliffs_delta = cd, direction = ifelse(med_kpc > med_ndm, 'KPC', 'NDM'))
  })
  mech_dt <- rbindlist(mech_rows)[order(p_value)]
  mech_dt[, fdr := p.adjust(p_value, 'fdr')]
  save_raw_table(mech_dt, 'raw_classic_final_biomarkers_kpc_ndm')

  cat('Saved RAW classic biomarker outputs to:', RAW_MODULAR_OUT, '\n')
}

main()
