source('scripts_raw/00_shared_utils_raw.R')

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

main <- function() {

  rep_ds <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_replicate_dataset.rds'))
  st_ds  <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_peak_strain_treatment_dataset.rds'))

  dists_rep <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_classic_distances_replicate.rds'))
  dists_st  <- readRDS(file.path(RAW_MODULAR_OBJ, 'raw_classic_distances_strain_treatment.rds'))

  datasets <- list(
    replicate        = list(meta = rep_ds$metadata, dists = dists_rep),
    strain_treatment = list(meta = st_ds$metadata,  dists = dists_st)
  )

  all_perm <- list()

  for (mode in names(datasets)) {
    meta  <- datasets[[mode]]$meta
    dists <- datasets[[mode]]$dists

    keep <- !is.na(meta$label_res_sen) &
            !is.na(meta$mechanism)    &
            !is.na(meta$source_lab)   &
            !is.na(meta$treatment)
    meta_c <- as.data.frame(meta[keep, ])
    cat(sprintf('[%s] Complete cases: %d / %d\n', mode, sum(keep), nrow(meta)))

    if (nrow(meta_c) < 10) {
      cat(sprintf('[%s] Too few complete cases, skipping.\n', mode))
      next
    }

    keep_idx <- which(keep)

    for (nm in names(dists)) {
      d   <- dists[[nm]]
      d_c <- as.dist(as.matrix(d)[keep_idx, keep_idx, drop = FALSE])

      terms <- c('label_res_sen', 'mechanism', 'source_lab', 'treatment')
      terms <- terms[vapply(terms, function(v)
        v %in% names(meta_c) && length(unique(na.omit(meta_c[[v]]))) > 1, logical(1))]
      frm <- as.formula(paste('d_c ~', paste(terms, collapse = ' + ')))

      result <- tryCatch(
        adonis2(frm, data = meta_c, permutations = 999),
        error = function(e) { cat(sprintf('  [%s %s] adonis2 error: %s\n', mode, nm, e$message)); NULL }
      )

      if (!is.null(result)) {
        all_perm[[length(all_perm) + 1]] <- data.table(
          group_mode     = mode,
          representation = nm,
          task           = 'multiclass',
          term           = rownames(as.data.frame(result)),
          as.data.table(as.data.frame(result), keep.rownames = FALSE)
        )
        cat(sprintf('  [%s %s] Model R2=%.4f F=%.2f p=%s\n',
                    mode, nm,
                    result['Model', 'R2'],
                    result['Model', 'F'],
                    result['Model', 'Pr(>F)']))
      }
    }
  }

  if (length(all_perm)) {
    out <- rbindlist(all_perm, fill = TRUE)
    save_raw_table(out, 'raw_classic_permanova_multiclass')
    cat('Saved: raw_classic_permanova_multiclass.csv\n')
    print(out[term == 'Model', .(group_mode, representation, task, R2, F, `Pr(>F)`)])
  } else {
    cat('No results to save.\n')
  }
}

main()
