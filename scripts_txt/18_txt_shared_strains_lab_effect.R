source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
})

aggregate_by_strain_treatment_lab <- function(mat, meta) {
  keys <- paste(meta$strain_id, meta$treatment, meta$source_lab, sep = "__")
  split_idx <- split(seq_len(nrow(meta)), keys)
  agg <- t(vapply(split_idx, function(idx) colMeans(mat[idx, , drop = FALSE]), numeric(ncol(mat))))
  rownames(agg) <- names(split_idx)

  det <- t(vapply(split_idx, function(idx) {
    as.numeric(colMeans((mat[idx, , drop = FALSE] > 0) * 1) >= 0.5)
  }, numeric(ncol(mat))))
  rownames(det) <- names(split_idx)
  colnames(det) <- colnames(mat)

  meta_agg <- rbindlist(lapply(names(split_idx), function(k) {
    idx <- split_idx[[k]]
    parts <- strsplit(k, "__")[[1]]
    data.table(
      sample     = k,
      strain_id  = parts[1],
      treatment  = parts[2],
      source_lab = parts[3],
      mechanism  = names(sort(table(meta$mechanism[idx]), decreasing = TRUE))[1],
      n_reps     = length(idx)
    )
  }))
  list(matrix = agg, detected = det, metadata = meta_agg)
}

cosine_dist_mat <- function(mat) {
  norms <- sqrt(rowSums(mat^2))
  norms[norms == 0] <- 1e-10
  mat_n <- mat / norms
  d <- 1 - tcrossprod(mat_n)
  d[d < 0] <- 0
  as.dist(d)
}

run_permanova_permadisp <- function(d, meta_a, rep_name, dist_name) {
  results <- list()

  perm <- adonis2(
    d ~ source_lab + treatment,
    data         = as.data.frame(meta_a),
    permutations = 999,
    strata       = meta_a$strain_id,
    by           = "margin"
  )
  dt <- as.data.table(perm, keep.rownames = "term")
  dt[, `:=`(representation = rep_name, distance = dist_name, test = "permanova_strata_strain")]
  results[[length(results) + 1]] <- dt

  disp_lab  <- betadisper(d, group = meta_a$source_lab)
  disp_perm <- permutest(disp_lab, permutations = 999)
  disp_dt   <- data.table(
    term           = rownames(disp_perm$tab),
    representation = rep_name,
    distance       = dist_name,
    test           = "permadisp_source_lab",
    Df             = disp_perm$tab$Df,
    SS             = disp_perm$tab$SS,
    MS             = disp_perm$tab$MS,
    F              = disp_perm$tab$F,
    `Pr(>F)`       = disp_perm$tab$`Pr(>F)`
  )
  results[[length(results) + 1]] <- disp_dt

  rbindlist(results, fill = TRUE)
}

main <- function() {
  meta_raw <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds"))
  mat_raw  <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))

  shared_ids <- meta_raw[, .(n_labs = uniqueN(source_lab)), by = strain_id][n_labs > 1]$strain_id
  cat("Cepas compartidas entre sitios:", length(shared_ids), "\n")

  keep    <- meta_raw$strain_id %in% shared_ids
  meta_sh <- meta_raw[keep]
  cat("Raw spectra:", nrow(meta_sh), "\n")

  agg <- aggregate_by_strain_treatment_lab(mat_raw[keep, , drop = FALSE], meta_sh)
  agg_meta <- agg$metadata
  mat_agg  <- agg$matrix
  cat("Aggregated samples:", nrow(agg_meta), "\n")
  print(agg_meta[, .N, by = .(source_lab, treatment)])
  save_table(agg_meta, "shared_strains_agg_metadata")

  agg_mats <- make_representations(mat_agg)
  agg_mats$binary <- agg$detected

  dist_config <- list(
    intensity_raw = "cosine",
    intensity_tic = "bray",
    intensity_max = "euclidean",
    intensity_clr = "euclidean",
    binary        = "jaccard"
  )

  results <- list()

  for (nm in names(dist_config)) {
    mat_a <- agg_mats[[nm]]
    for (dist_name in dist_config[[nm]]) {
      cat("Running:", nm, "/", dist_name, "\n")
      d <- switch(dist_name,
        cosine     = cosine_dist_mat(mat_a),
        bray       = vegdist(mat_a, method = "bray"),
        euclidean  = dist(mat_a, method = "euclidean"),
        jaccard    = vegdist(mat_a, method = "jaccard", binary = TRUE)
      )
      results[[length(results) + 1]] <- run_permanova_permadisp(d, agg_meta, nm, dist_name)

      pcoa <- cmdscale(d, k = 2, eig = TRUE)
      eigs <- pcoa$eig; eigs[eigs < 0] <- 0
      pct  <- round(eigs[1:2] / sum(eigs) * 100, 1)
      pcoa_dt <- data.table(
        PC1        = pcoa$points[, 1],
        PC2        = pcoa$points[, 2],
        source_lab = agg_meta$source_lab,
        treatment  = agg_meta$treatment,
        strain_id  = agg_meta$strain_id,
        mechanism  = agg_meta$mechanism
      )
      p <- ggplot(pcoa_dt, aes(PC1, PC2, color = source_lab, shape = treatment)) +
        geom_point(size = 2.5, alpha = 0.85) +
        geom_line(aes(group = paste(strain_id, treatment)), color = "grey65", linewidth = 0.3) +
        scale_color_manual(values = c(ALVAREZ = "#E64B35", MALBRAN = "#4DBBD5")) +
        theme_bw() +
        labs(
          title = paste0("PCoA (", dist_name, ") — ", nm, " — ",
                                 uniqueN(agg_meta$strain_id), " shared strains"),
          x     = paste0("PC1 (", pct[1], "%)"),
          y     = paste0("PC2 (", pct[2], "%)")
        )
      save_plot(p, paste0("shared_strains_pcoa_", nm, "_", dist_name), width = 7, height = 5)
    }
  }

  out <- rbindlist(results, fill = TRUE)
  save_table(out, "shared_strains_lab_effect_permanova")

  cat("\n=== Summary: source_lab R² and p across representations ===\n")
  print(out[test == "permanova_strata_strain" & term == "source_lab",
            .(representation, distance, R2, F, `Pr(>F)`)])
  cat("\nSaved outputs to:", TXT_MODULAR_TAB, "\n")
}

main()
