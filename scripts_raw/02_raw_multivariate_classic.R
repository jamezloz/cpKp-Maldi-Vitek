source('scripts_raw/00_shared_utils_raw.R')

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

main <- function() {
  need_pkg('vegan')
  rep_ds <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_replicate_dataset.rds'))
  strain_ds <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_strain_mean_dataset.rds'))
  st_ds <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_strain_treatment_dataset.rds'))

  datasets <- list(replicate = rep_ds, strain_mean = strain_ds, strain_treatment = st_ds)
  all_perm <- list()
  all_disp <- list()

  for (mode in names(datasets)) {
    ds <- datasets[[mode]]
    meta <- ds$metadata
    reps <- make_raw_peak_representations(ds$matrix, ds$detected)

    dists <- list(
      peak_intensity_raw_cosine    = cosine_dist(reps$peak_intensity_raw),
      peak_intensity_tic_bray      = vegan::vegdist(reps$peak_intensity_tic, method = 'bray'),
      peak_intensity_max_euclidean = dist(reps$peak_intensity_max),
      peak_intensity_clr_euclidean = dist(reps$peak_intensity_clr),
      peak_binary_jaccard          = jaccard_dist(reps$peak_binary)
    )

    save_raw_object(dists, paste0('raw_classic_distances_', mode))

    adjust_treatment <- mode %in% c('replicate', 'strain_treatment')
    for (nm in names(dists)) {
      d <- dists[[nm]]
      if (sum(!is.na(meta$label_res_sen)) > 3) {
        frm_res <- if (adjust_treatment) d ~ label_res_sen + treatment else d ~ label_res_sen
        perm_res <- vegan::adonis2(frm_res, data = as.data.frame(meta), permutations = 999, by = 'margin')
        all_perm[[length(all_perm) + 1]] <- data.table(group_mode = mode, representation = nm, task = 'res_sen', term = rownames(as.data.frame(perm_res)), as.data.table(perm_res, keep.rownames = FALSE))
      }
      keep_mech <- meta$mechanism %in% c('KPC', 'NDM')
      if (sum(keep_mech) > 3) {
        d_mech <- as.dist(as.matrix(d)[keep_mech, keep_mech, drop = FALSE])
        meta_mech <- as.data.frame(meta[keep_mech])
        frm_mech <- if (adjust_treatment) d_mech ~ mechanism + treatment else d_mech ~ mechanism
        perm_mech <- vegan::adonis2(frm_mech, data = meta_mech, permutations = 999, by = 'margin')
        all_perm[[length(all_perm) + 1]] <- data.table(group_mode = mode, representation = nm, task = 'kpc_ndm', term = rownames(as.data.frame(perm_mech)), as.data.table(perm_mech, keep.rownames = FALSE))
      }
      bd <- tryCatch(vegan::betadisper(d, meta$label_res_sen), error = function(e) NULL)
      if (!is.null(bd)) {
        all_disp[[length(all_disp) + 1]] <- data.table(group_mode = mode, representation = nm, task = 'res_sen', F_value = anova(bd)$`F value`[1], p_value = anova(bd)$`Pr(>F)`[1])
      }

      pcoa <- cmdscale(d, k = 2, eig = TRUE)
      df <- data.frame(Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2], label_res_sen = meta$label_res_sen, mechanism = meta$mechanism)
      p <- ggplot(df, aes(Axis1, Axis2, color = label_res_sen, shape = mechanism)) + geom_point(alpha = 0.8) + theme_bw() + labs(title = paste('RAW classic', mode, nm))
      save_raw_plot(p, paste0('raw_', mode, '_', nm, '_pcoa'), width = 8, height = 6)
    }
  }

  if (length(all_perm)) save_raw_table(rbindlist(all_perm, fill = TRUE), 'raw_classic_permanova')
  if (length(all_disp)) save_raw_table(rbindlist(all_disp, fill = TRUE), 'raw_classic_permadisp')
  cat('Saved RAW classic multivariate outputs to:', RAW_MODULAR_OUT, '\n')
}

main()
