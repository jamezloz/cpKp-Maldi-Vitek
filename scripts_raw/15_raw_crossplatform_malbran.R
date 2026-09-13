source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
})

build_pooled <- function() {
  meta <- build_raw_metadata('pooled')
  cat("Espectros:", nrow(meta), " | ")
  print(table(meta$source_lab, meta$label_res_sen))
  specs <- lapply(meta$path, preprocess_raw_spectrum, baseline = 'tophat', normalize = FALSE)
  names(specs) <- meta$sample
  peak_obj <- build_peak_matrix_raw(specs)
  list(matrix = peak_obj$matrix[meta$sample, , drop = FALSE], metadata = meta)
}

crossplatform_distances <- function(pooled) {
  agg <- aggregate_raw_by_group(pooled$matrix, pooled$metadata, 'strain_treatment')
  m <- as.data.table(agg$metadata)
  X <- agg$matrix
  D <- as.matrix(cosine_dist(X))

  shared <- intersect(m[source_lab == "ALVAREZ", unique(strain_id)],
                      m[source_lab == "MALBRAN", unique(strain_id)])
  cat("\nCepas SEN medidas en ambos equipos:", length(shared), "\n")

  rows <- list()
  for (s in shared) {
    ia <- which(m$strain_id == s & m$source_lab == "ALVAREZ")
    im <- which(m$strain_id == s & m$source_lab == "MALBRAN")
    if (!length(ia) || !length(im)) next
    cross_same <- mean(D[ia, im])
    others_a <- which(m$strain_id != s & m$source_lab == "ALVAREZ" & m$label_res_sen == "SEN")
    others_m <- which(m$strain_id != s & m$source_lab == "MALBRAN")
    rows[[length(rows) + 1]] <- data.table(
      strain_id = s,
      intra_cepa_cross_sitio = cross_same,
      inter_cepa_alvarez = if (length(others_a)) mean(D[ia, others_a]) else NA_real_,
      inter_cepa_malbran = if (length(others_m)) mean(D[im, others_m]) else NA_real_)
    }
  res <- rbindlist(rows, fill = TRUE)
  save_raw_table(res, "raw_crossplatform_distances")

  summ <- data.table(
    n_cepas = nrow(res),
    d_misma_cepa_otro_equipo = round(mean(res$intra_cepa_cross_sitio, na.rm = TRUE), 4),
    d_otras_cepas_alvarez = round(mean(res$inter_cepa_alvarez, na.rm = TRUE), 4),
    d_otras_cepas_malbran = round(mean(res$inter_cepa_malbran, na.rm = TRUE), 4))

  wt <- wilcox.test(res$intra_cepa_cross_sitio, res$inter_cepa_alvarez, paired = TRUE)
  summ[, wilcoxon_p := wt$p.value]
  save_raw_table(summ, "raw_crossplatform_summary")
  cat("\nReproducibilidad inter-plataforma:\n"); print(summ)
  invisible(res)
}

external_specificity <- function(pooled) {
  agg <- aggregate_raw_by_group(pooled$matrix, pooled$metadata, 'strain_treatment')
  m <- as.data.table(agg$metadata)
  X <- agg$matrix

  tr <- which(m$source_lab == "ALVAREZ" & !is.na(m$label_res_sen))
  te <- which(m$source_lab == "MALBRAN" & !is.na(m$label_res_sen))
  shared <- intersect(m$strain_id[tr], m$strain_id[te])

  tr_disj <- setdiff(tr, which(m$strain_id %in% shared))
  disj_viable <- uniqueN(m$label_res_sen[tr_disj]) >= 2
  cat("Cepas SEN de Alvarez que NO estan en Malbran:",
      uniqueN(m$strain_id[tr_disj][m$label_res_sen[tr_disj] == "SEN"]),
      "-> variante disjunta", if (disj_viable) "viable" else "imposible", "\n")

  run <- function(train_idx, label) {
    fit <- randomForest(x = X[train_idx, , drop = FALSE],
                        y = factor(m$label_res_sen[train_idx]), ntree = 500)
    pred <- as.character(predict(fit, X[te, , drop = FALSE]))
    obs <- m$label_res_sen[te]
    data.table(analisis = label,
               cepas_train = uniqueN(m$strain_id[train_idx]),
               muestras_test = length(te),
               cepas_test = uniqueN(m$strain_id[te]),
               especificidad = round(mean(pred[obs == "SEN"] == "SEN"), 3),
               n_falsos_positivos = sum(pred[obs == "SEN"] != "SEN"))
  }

  res <- run(tr, "alvarez_completo")
  if (disj_viable) {
    res <- rbindlist(list(res, run(tr_disj, "alvarez_sin_cepas_compartidas")))
  }
  save_raw_table(res, "raw_external_specificity_malbran")
  cat("\nValidacion externa en Malbran (todas sensibles -> solo especificidad):\n")
  print(res)
  invisible(res)
}

site_effect_within_sen <- function(pooled) {
  safe <- requireNamespace("vegan", quietly = TRUE)
  if (!safe) return(NULL)
  agg <- aggregate_raw_by_group(pooled$matrix, pooled$metadata, 'strain_treatment')
  m <- as.data.table(agg$metadata)
  sen <- which(m$label_res_sen == "SEN")
  X <- agg$matrix[sen, , drop = FALSE]
  ms <- droplevels(as.data.frame(m[sen]))
  d <- cosine_dist(X)

  fit <- vegan::adonis2(d ~ source_lab + treatment, data = ms,
                        permutations = 999, by = "margin")
  disp <- vegan::betadisper(d, ms$source_lab)
  res <- data.table(term = rownames(fit), as.data.table(fit, keep.rownames = NULL))
  res[, permadisp_F := NA_real_][, permadisp_p := NA_real_]
  res[term == "source_lab", `:=`(permadisp_F = anova(disp)$`F value`[1],
                                 permadisp_p = anova(disp)$`Pr(>F)`[1])]
  save_raw_table(res, "raw_site_effect_within_sen")
  cat("\nEfecto de sitio con fenotipo constante (solo sensibles):\n")
  print(res[term %in% c("source_lab", "treatment"),
            .(term, Df, R2 = round(R2, 4), p = `Pr(>F)`)])
  cat("n muestras:", length(sen), " cepas:", uniqueN(ms$strain_id), "\n")
  invisible(res)
}

main <- function() {
  set.seed(42)
  pooled <- build_pooled()
  saveRDS(pooled, file.path(RAW_MODULAR_OBJ, "raw_pooled_dataset.rds"))
  crossplatform_distances(pooled)
  site_effect_within_sen(pooled)
  external_specificity(pooled)
  cat("\nSalidas en:", RAW_MODULAR_TAB, "\n")
}

main()
