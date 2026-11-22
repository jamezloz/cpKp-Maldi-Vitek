source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

main <- function() {
  meta_raw <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds"))
  mat_raw <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))

  # Strains measured at both sites
  shared_strains <- meta_raw[, .(n_labs = uniqueN(source_lab)), by = strain_id][n_labs > 1]$strain_id
  # Index by position: 189 peak-list files share their name across sites
  stopifnot(identical(rownames(mat_raw), as.character(meta_raw$sample)))
  ix <- which(meta_raw$strain_id %in% shared_strains)
  shared <- meta_raw[ix]
  save_table(shared, "shared_strains_metadata_raw")

  # Build pairwise within-strain / cross-lab distances on raw spectra
  d <- as.matrix(cosine_dist(mat_raw[ix, , drop = FALSE]))
  rows <- list()
  for (i in seq_len(nrow(shared))) {
    for (j in seq_len(nrow(shared))) {
      if (j <= i) next
      rows[[length(rows) + 1]] <- data.table(
        sample_a = shared$sample[i],
        sample_b = shared$sample[j],
        strain_a = shared$strain_id[i],
        strain_b = shared$strain_id[j],
        same_strain = shared$strain_id[i] == shared$strain_id[j],
        same_treatment = shared$treatment[i] == shared$treatment[j],
        lab_a = shared$source_lab[i],
        lab_b = shared$source_lab[j],
        cross_lab = shared$source_lab[i] != shared$source_lab[j],
        cosine_distance = d[i, j]
      )
    }
  }
  pair_dt <- rbindlist(rows)
  save_table(pair_dt, "shared_strains_pairwise_distances")

  # Summary
  summary_dt <- pair_dt[, .(
    mean_cosine_distance = mean(cosine_distance, na.rm = TRUE),
    median_cosine_distance = median(cosine_distance, na.rm = TRUE),
    n = .N
  ), by = .(same_strain, same_treatment, cross_lab)]
  save_table(summary_dt, "shared_strains_distance_summary")

  if (nrow(pair_dt)) {
    p <- ggplot(pair_dt, aes(x = interaction(same_strain, cross_lab), y = cosine_distance, fill = same_strain)) +
      geom_boxplot(outlier.alpha = 0.2) +
      theme_bw() +
      labs(x = "same_strain.cross_lab", y = "Cosine distance", title = "Shared strains across MALBRAN and ALVAREZ")
    save_plot(p, "shared_strains_cosine_distance_boxplot", width = 8, height = 5)
  }

  cat("Saved shared-strain outputs to:", TXT_MODULAR_OUT, "\n")
}

main()
