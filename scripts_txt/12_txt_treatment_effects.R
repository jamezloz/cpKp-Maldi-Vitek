source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

extract_adonis_term <- function(ad, representation, model_label, term_name) {
  dt <- as.data.table(ad, keep.rownames = "term")
  dt <- dt[term == term_name]
  if (nrow(dt) == 0) return(NULL)
  data.table(
    representation = representation,
    model = model_label,
    term = term_name,
    df = dt$Df,
    sumsq = dt$SumOfSqs,
    r2 = dt$R2,
    f_stat = dt$F,
    p_value = dt$`Pr(>F)`
  )
}

main <- function() {
  safe_require("vegan")

  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))

  canonical_reps <- c("intensity_raw_cosine", "intensity_tic_bray",
                      "intensity_max_euclidean", "intensity_clr_euclidean",
                      "binary_jaccard")
  dist_files <- file.path(TXT_MODULAR_OBJ, paste0("distance_", canonical_reps, ".rds"))
  stopifnot(all(file.exists(dist_files)))
  dists <- setNames(lapply(dist_files, readRDS), canonical_reps)

  meta_treat <- copy(meta[!is.na(treatment) & !is.na(label_res_sen) & !is.na(source_lab)])
  keep_idx <- meta$sample %in% meta_treat$sample

  rows <- list()
  for (nm in names(dists)) {
    d_full <- as.matrix(dists[[nm]])
    d <- as.dist(d_full[keep_idx, keep_idx, drop = FALSE])

    ad_full <- vegan::adonis2(
      d ~ treatment + label_res_sen + source_lab,
      data = as.data.frame(meta_treat),
      permutations = 999,
      by = "margin"
    )
    rows[[length(rows) + 1]] <- extract_adonis_term(ad_full, nm, "treatment+label_res_sen+source_lab", "treatment")
    rows[[length(rows) + 1]] <- extract_adonis_term(ad_full, nm, "treatment+label_res_sen+source_lab", "source_lab")
    rows[[length(rows) + 1]] <- extract_adonis_term(ad_full, nm, "treatment+label_res_sen+source_lab", "label_res_sen")

    ad_strata <- vegan::adonis2(
      d ~ label_res_sen + source_lab,
      data = as.data.frame(meta_treat),
      permutations = 999,
      strata = meta_treat$treatment,
      by = "margin"
    )
    rows[[length(rows) + 1]] <- extract_adonis_term(ad_strata, nm, "label_res_sen+source_lab|strata=treatment", "label_res_sen")
    rows[[length(rows) + 1]] <- extract_adonis_term(ad_strata, nm, "label_res_sen+source_lab|strata=treatment", "source_lab")
  }
  out <- rbindlist(rows, fill = TRUE)
  save_table(out, "treatment_platform_permanova")

  raw_meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds"))
  raw_mat <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))
  raw_cos <- as.matrix(cosine_dist(raw_mat))

  pair_rows <- list()
  for (sid in unique(raw_meta$strain_id)) {
    g <- raw_meta[strain_id == sid & !is.na(treatment)]
    if (nrow(g) == 0 || !all(c("COL", "EXT") %in% unique(g$treatment))) next

    idx_col <- which(raw_meta$sample %in% g[treatment == "COL"]$sample)
    idx_ext <- which(raw_meta$sample %in% g[treatment == "EXT"]$sample)
    vals <- as.vector(raw_cos[idx_col, idx_ext, drop = FALSE])
    mech <- names(sort(table(g$mechanism), decreasing = TRUE))[1]
    labs <- unique(g$source_lab)
    lab_label <- if (length(labs) == 1) labs else "MIXED"

    pair_rows[[length(pair_rows) + 1]] <- data.table(
      strain_id = sid,
      n_col = length(idx_col),
      n_ext = length(idx_ext),
      mean_col_ext_distance = mean(vals, na.rm = TRUE),
      median_col_ext_distance = median(vals, na.rm = TRUE),
      mechanism = mech,
      source_lab = lab_label,
      label_res_sen = names(sort(table(g$label_res_sen), decreasing = TRUE))[1]
    )
  }

  if (length(pair_rows)) {
    pair_dt <- rbindlist(pair_rows, fill = TRUE)
    save_table(pair_dt, "treatment_col_ext_pair_summary")

    p1 <- ggplot(pair_dt, aes(mechanism, mean_col_ext_distance, fill = source_lab)) +
      geom_boxplot(outlier.alpha = 0.2) +
      theme_bw() +
      labs(title = "Within-strain COL vs EXT distance", y = "Mean cosine distance", x = "Mechanism")
    save_plot(p1, "treatment_col_ext_distance_boxplot", width = 8, height = 5)

    p2 <- ggplot(pair_dt, aes(label_res_sen, mean_col_ext_distance, fill = source_lab)) +
      geom_boxplot(outlier.alpha = 0.2) +
      theme_bw() +
      labs(title = "Within-strain COL vs EXT distance by phenotype", y = "Mean cosine distance", x = "Phenotype")
    save_plot(p2, "treatment_col_ext_distance_by_phenotype", width = 7, height = 5)
  }

  cat("Saved treatment/platform outputs to:", TXT_MODULAR_OUT, "\n")
}

main()
