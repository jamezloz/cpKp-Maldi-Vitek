source("scripts_txt/00_shared_utils.R")

cliffs_delta_manual <- function(x, y) {
  cmp <- outer(x, y, FUN = function(a, b) sign(a - b))
  sum(cmp) / (length(x) * length(y))
}

stability_selection_proxy <- function(mat, y, strain, n_boot = 100) {

  freq <- setNames(numeric(ncol(mat)), colnames(mat))
  strains <- unique(strain)
  for (b in seq_len(n_boot)) {
    boot_strains <- sample(strains, replace = TRUE)
    idx <- unlist(lapply(boot_strains, function(s) which(strain == s)))
    boot_mat <- mat[idx, , drop = FALSE]
    boot_y <- y[idx]
    if (length(unique(boot_y)) < 2) next
    pvals <- apply(boot_mat, 2, function(v)
      tryCatch(wilcox.test(v ~ boot_y, exact = FALSE)$p.value, error = function(e) 1))
    fdr <- p.adjust(pvals, method = "fdr")

    sig <- which(fdr < 0.05)
    freq[sig] <- freq[sig] + 1
  }
  data.table(feature = names(freq), selection_freq = freq / n_boot)[order(-selection_freq)]
}

collapse_lab_replicates <- function(mat, meta) {
  key <- paste(meta$strain_id, meta$treatment, sep = "__")
  idx <- split(seq_len(nrow(meta)), key)
  m <- t(vapply(idx, function(i) colMeans(mat[i, , drop = FALSE]), numeric(ncol(mat))))
  md <- rbindlist(lapply(idx, function(i) {
    g <- meta[i]
    data.table(sample = paste(g$strain_id[1], g$treatment[1], sep = "__"),
               strain_id = g$strain_id[1], treatment = g$treatment[1],
               mechanism = g$mechanism[1], mechanism_raw = g$mechanism_raw[1],
               label_res_sen = g$label_res_sen[1])
  }))
  rownames(m) <- md$sample
  list(matrix = m, metadata = md)
}

univariate_tests <- function(mat_f, meta, level) {
  kr_p <- apply(mat_f, 2, function(v) kruskal.test(v ~ factor(meta$mechanism_raw))$p.value)
  kr <- data.table(level = level, feature = colnames(mat_f), p_value = kr_p,
                   fdr = p.adjust(kr_p, method = "fdr"))[order(fdr, p_value)]

  contrasts <- list(c("SEN", "KPC"), c("SEN", "NDM"), c("KPC", "NDM"), c("SEN", "KPC-NDM"))
  cliff_rows <- list()
  for (pair in contrasts) {
    idx <- meta$mechanism_raw %in% pair
    subm <- mat_f[idx, , drop = FALSE]
    suby <- factor(meta$mechanism_raw[idx], levels = pair)
    if (length(unique(suby)) < 2) next
    cliff_rows[[length(cliff_rows) + 1]] <- rbindlist(lapply(seq_len(ncol(subm)), function(j) {
      data.table(level = level, comparison = paste(pair, collapse = "_vs_"),
                 feature = colnames(subm)[j],
                 cliffs_delta = cliffs_delta_manual(subm[suby == pair[1], j], subm[suby == pair[2], j]))
    }))
  }
  list(kruskal = kr, cliff = rbindlist(cliff_rows, fill = TRUE))
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  mat <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_aggregated.rds"))

  bin <- mat > 0
  res_i <- meta$label_res_sen == "RES"
  prev_max <- pmax(colMeans(bin[res_i, , drop = FALSE]),
                   colMeans(bin[!res_i, , drop = FALSE]))

  keep <- prev_max >= 0.05
  mat_f <- mat[, keep, drop = FALSE]
  cat(sprintf("  Features: %d -> %d (prevalencia maxima por clase >= 5%%)\n",
              ncol(mat), ncol(mat_f)))

  prot_ds <- collapse_lab_replicates(mat_f, as.data.table(meta))
  prim <- univariate_tests(prot_ds$matrix, prot_ds$metadata, "cepa_x_protocolo")
  sens <- univariate_tests(mat_f, as.data.table(meta), "cepa_x_protocolo_x_sitio")

  save_table(rbindlist(list(prim$kruskal, sens$kruskal)), "biomarkers_global_kruskal")
  save_table(rbindlist(list(prim$cliff, sens$cliff)), "biomarkers_global_cliff")

  cat("Features con FDR < 0.05 (Kruskal):\n")
  cat("  cepa x protocolo (primario) (n =", nrow(prot_ds$matrix), "):",
      prim$kruskal[fdr < 0.05, .N], "\n")
  cat("  cepa x protocolo x sitio    (n =", nrow(mat_f), "):",
      sens$kruskal[fdr < 0.05, .N], "\n")

  pm <- prot_ds$metadata
  stab_res <- stability_selection_proxy(prot_ds$matrix, factor(pm$label_res_sen), pm$strain_id)
  save_table(stab_res, "stability_selection_res_sen")

  keep_mech <- pm$mechanism %in% c("KPC", "NDM")
  stab_mech <- stability_selection_proxy(prot_ds$matrix[keep_mech, , drop = FALSE],
                                         factor(pm$mechanism[keep_mech]), pm$strain_id[keep_mech])
  save_table(stab_mech, "stability_selection_kpc_ndm")

  cat("Saved global biomarker outputs to:", TXT_MODULAR_OUT, "\n")
}

main()
