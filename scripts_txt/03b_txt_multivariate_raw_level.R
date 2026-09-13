source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(vegan)
  library(cluster)
})

canonical_dist <- function(rep_name, mat) {
  switch(rep_name,
    intensity_raw = cosine_dist(mat),
    intensity_tic = vegan::vegdist(mat, method = "bray"),
    intensity_max = dist(mat, method = "euclidean"),
    intensity_clr = dist(mat, method = "euclidean"),
    binary        = vegan::vegdist(mat, method = "jaccard", binary = TRUE)
  )
}

run_permanova_permadisp <- function(rep_name, d, meta) {
  out <- list()

  perm <- adonis2(
    d ~ label_res_sen + mechanism + source_lab + treatment,
    data = meta, permutations = 999, by = "terms"
  )
  dt <- as.data.table(perm, keep.rownames = "term")
  dt[, `:=`(representation = rep_name, level = "raw_spectrum", test = "permanova")]
  out[[1]] <- dt

  perm_strata <- adonis2(
    d ~ label_res_sen + mechanism + source_lab,
    data = meta, permutations = 999, strata = meta$treatment, by = "terms"
  )
  dt2 <- as.data.table(perm_strata, keep.rownames = "term")
  dt2[, `:=`(representation = rep_name, level = "raw_spectrum", test = "permanova_stratified_treatment")]
  out[[2]] <- dt2

  perm_res <- adonis2(d ~ label_res_sen + source_lab + treatment,
                      data = meta, permutations = 999, by = "margin")
  dt3 <- as.data.table(perm_res, keep.rownames = "term")
  dt3[, `:=`(representation = rep_name, level = "raw_spectrum", test = "margin_res_sen")]
  out[[3]] <- dt3

  keep_kn <- meta$mechanism %in% c("KPC", "NDM")
  if (sum(keep_kn) > 20 && uniqueN(meta$mechanism[keep_kn]) > 1) {
    d_mat <- as.matrix(d)
    d_kn  <- as.dist(d_mat[keep_kn, keep_kn, drop = FALSE])
    meta_kn <- droplevels(meta[keep_kn, ])
    fit <- try(adonis2(d_kn ~ mechanism + source_lab + treatment,
                       data = meta_kn, permutations = 999, by = "margin"),
              silent = TRUE)
    if (!inherits(fit, "try-error")) {
      dt4 <- as.data.table(fit, keep.rownames = "term")
      dt4[, `:=`(representation = rep_name, level = "raw_spectrum", test = "margin_kpc_ndm")]
      out[[4]] <- dt4
    }
  }

  disp <- betadisper(d, group = meta$label_res_sen)
  disp_dt <- as.data.table(anova(disp), keep.rownames = "term")
  disp_dt[, `:=`(representation = rep_name, level = "raw_spectrum", test = "permadisp_label_res_sen")]
  out[[length(out) + 1]] <- disp_dt

  sil <- cluster::silhouette(as.integer(as.factor(meta$label_res_sen)), d)
  sil_dt <- as.data.table(sil[, c("cluster", "neighbor", "sil_width")])
  sil_dt[, `:=`(representation = rep_name, level = "raw_spectrum")]
  save_table(sil_dt, paste0("raw_level_silhouette_", rep_name))

  rbindlist(out, fill = TRUE)
}

main <- function() {
  meta_raw <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds"))
  mat_raw  <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))
  meta_df  <- as.data.frame(meta_raw)

  cat("Raw spectra:", nrow(meta_raw), "\n")

  reps_all <- make_representations(mat_raw)
  rep_list <- list(
    intensity_raw = reps_all$intensity_raw,
    intensity_tic = reps_all$intensity_tic,
    intensity_max = reps_all$intensity_max,
    intensity_clr = reps_all$intensity_clr,
    binary        = reps_all$binary
  )

  results <- list()
  for (nm in names(rep_list)) {
    cat("Computing distance:", nm, "\n")
    d <- canonical_dist(nm, rep_list[[nm]])
    cat("Running PERMANOVA:", nm, "\n")
    results[[nm]] <- run_permanova_permadisp(nm, d, meta_df)
  }

  out <- rbindlist(results, fill = TRUE)
  save_table(out, "permanova_summary_raw_level")

  raw_terms <- out[test == "permanova" & !term %in% c("Residual", "Total"),
                   .(representation, term, R2, F, `Pr(>F)`)]
  save_table(raw_terms, "permanova_raw_vs_aggregated_comparison")

  model_r2 <- function(dt) {
    m <- dt[test == "permanova" & term == "Model", .(representation, R2)]
    if (nrow(m)) return(m)
    dt[test == "permanova" & term == "Residual", .(representation, R2 = 1 - R2)]
  }
  agg_perm <- fread(file.path(TXT_MODULAR_TAB, "permanova_summary.csv"))
  agg_model <- model_r2(agg_perm); setnames(agg_model, "R2", "R2_aggregated")
  raw_model <- model_r2(out);      setnames(raw_model, "R2", "R2_raw")
  comp <- merge(agg_model, raw_model, by = "representation", all.y = TRUE)
  stopifnot(nrow(comp) > 0)
  save_table(comp, "permanova_model_r2_comparison")

  cat("\n=== Marginal res_sen (ajustado por sitio+protocolo), nivel replica ===\n")
  print(out[test == "margin_res_sen" & term == "label_res_sen",
            .(representation, R2, F, `Pr(>F)`)])
  cat("\n=== Marginal KPC/NDM (ajustado por sitio+protocolo), nivel replica ===\n")
  print(out[test == "margin_kpc_ndm" & term == "mechanism",
            .(representation, R2, F, `Pr(>F)`)])
  cat("\n=== R2 del modelo completo (composite): agregado vs replica ===\n")
  print(comp)
  cat("\nSaved outputs to:", TXT_MODULAR_TAB, "\n")
}

main()
