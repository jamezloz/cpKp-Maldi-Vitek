source("scripts_txt/00_shared_utils.R")

safe_adonis <- function(d, meta, formula_text, strata = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NULL)
  frm <- as.formula(formula_text)
  vegan::adonis2(d ~ label_res_sen + mechanism + source_lab + treatment, data = meta, permutations = 999, strata = strata)
}

run_permanova_suite <- function(name, d, meta) {
  out <- list()
  if (requireNamespace("vegan", quietly = TRUE)) {

    perm <- vegan::adonis2(d ~ label_res_sen + mechanism + source_lab + treatment, data = meta, permutations = 999)
    out[[length(out) + 1]] <- data.table(representation = name, test = "permanova", term = rownames(perm), as.data.table(perm, keep.rownames = NULL))

    perm_res <- vegan::adonis2(d ~ label_res_sen + source_lab + treatment,
                               data = meta, permutations = 999, by = "margin")
    out[[length(out) + 1]] <- data.table(representation = name, test = "permanova_margin_res_sen",
                                         term = rownames(perm_res),
                                         as.data.table(perm_res, keep.rownames = NULL))

    res_only <- which(meta$label_res_sen == "RES")
    if (length(res_only) > 20 && uniqueN(meta$mechanism[res_only]) > 1) {
      d_res <- as.dist(as.matrix(d)[res_only, res_only, drop = FALSE])

      meta_res <- droplevels(as.data.frame(meta[res_only, ]))
      perm_mech <- try(vegan::adonis2(d_res ~ mechanism + source_lab + treatment,
                                      data = meta_res, permutations = 999, by = "margin"),
                       silent = TRUE)
      if (!inherits(perm_mech, "try-error")) {
        out[[length(out) + 1]] <- data.table(representation = name, test = "permanova_margin_mechanism",
                                             term = rownames(perm_mech),
                                             as.data.table(perm_mech, keep.rownames = NULL))
      }
    }
    perm_strata <- vegan::adonis2(d ~ label_res_sen + mechanism + source_lab, data = meta, permutations = 999, strata = meta$treatment)
    out[[length(out) + 1]] <- data.table(representation = name, test = "permanova_stratified_treatment", term = rownames(perm_strata), as.data.table(perm_strata, keep.rownames = NULL))
    disp <- vegan::betadisper(d, group = meta$label_res_sen)
    save_table(as.data.table(anova(disp), keep.rownames = "term"), paste0("permadisp_", name, "_anova"))
    save_object(disp, paste0(name, "_permadisp"))

    nmds <- vegan::metaMDS(d, k = 2, try = 20, trymax = 50, trace = FALSE)
    save_object(nmds, paste0(name, "_nmds"))

    pcoa <- ape::pcoa(d)
    save_object(pcoa, paste0(name, "_pcoa"))

    cap <- vegan::capscale(d ~ label_res_sen + source_lab + treatment, data = meta)
    save_object(cap, paste0(name, "_cap"))

    dbrda <- vegan::dbrda(d ~ label_res_sen + source_lab + treatment, data = meta)
    save_object(dbrda, paste0(name, "_dbrda"))
  }
  if (length(out)) rbindlist(out, fill = TRUE) else NULL
}

run_procrustes <- function(d1, d2) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NULL)
  ord1 <- cmdscale(d1, k = 3, eig = TRUE)
  ord2 <- cmdscale(d2, k = 3, eig = TRUE)
  fit <- vegan::procrustes(ord1$points, ord2$points)
  list(fit = fit, protest = vegan::protest(ord1$points, ord2$points, permutations = 999))
}

run_mantel <- function(d1, d2) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NULL)
  vegan::mantel(d1, d2, permutations = 999)
}

run_silhouette_proxy <- function(d, groups) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(NULL)
  grp <- as.integer(as.factor(groups))
  sil <- cluster::silhouette(grp, d)
  data.table(cluster = sil[, "cluster"], neighbor = sil[, "neighbor"], sil_width = sil[, "sil_width"])
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))

  canonical_reps <- c("intensity_raw_cosine", "intensity_tic_bray",
                      "intensity_max_euclidean", "intensity_clr_euclidean",
                      "binary_jaccard")
  dist_files <- file.path(TXT_MODULAR_OBJ, paste0("distance_", canonical_reps, ".rds"))
  stopifnot(all(file.exists(dist_files)))
  dists <- setNames(lapply(dist_files, readRDS), canonical_reps)

  perm_rows <- list()
  for (nm in names(dists)) {
    d <- dists[[nm]]
    res <- run_permanova_suite(nm, d, as.data.frame(meta))
    if (!is.null(res)) perm_rows[[length(perm_rows) + 1]] <- res

    pca <- prcomp(reps$intensity_raw, center = TRUE, scale. = TRUE)
    save_object(pca, paste0(nm, "_pca"))

    if (requireNamespace("ape", quietly = TRUE)) {
      nj_tree <- ape::nj(d)
      save_object(nj_tree, paste0(nm, "_nj"))
    }
    upgma <- hclust(d, method = "average")
    save_object(upgma, paste0(nm, "_upgma"))
    sil <- run_silhouette_proxy(d, meta$label_res_sen)
    if (!is.null(sil)) save_table(sil, paste0(nm, "_silhouette"))
  }

  if (length(names(dists)) >= 2) {
    combos <- combn(names(dists), 2, simplify = FALSE)
    pro_rows <- list()
    man_rows <- list()
    for (cb in combos) {
      pro <- run_procrustes(dists[[cb[1]]], dists[[cb[2]]])
      if (!is.null(pro)) {
        pro_rows[[length(pro_rows) + 1]] <- data.table(
          dist_a = cb[1],
          dist_b = cb[2],
          ss = pro$fit$ss,
          correlation = pro$protest$t0,
          p_value = pro$protest$signif
        )
      }
      man <- run_mantel(dists[[cb[1]]], dists[[cb[2]]])
      if (!is.null(man)) {
        man_rows[[length(man_rows) + 1]] <- data.table(
          dist_a = cb[1],
          dist_b = cb[2],
          statistic = man$statistic,
          p_value = man$signif
        )
      }
    }
    if (length(pro_rows)) save_table(rbindlist(pro_rows), "procrustes_summary")
    if (length(man_rows)) save_table(rbindlist(man_rows), "mantel_summary")
  }

  if (length(perm_rows)) save_table(rbindlist(perm_rows, fill = TRUE), "permanova_summary")
  cat("Saved multivariate outputs to:", TXT_MODULAR_OUT, "\n")
}

main()
