source("scripts_figures/00_fig_utils.R")
source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(vegan)
})

TXT_OBJ <- TXT_MODULAR_OBJ

cosine_dist_mat <- function(mat) {
  norms <- sqrt(rowSums(mat^2)); norms[norms == 0] <- 1e-10
  mat_n <- mat / norms
  d <- 1 - tcrossprod(mat_n); d[d < 0] <- 0
  as.dist(d)
}

meta_raw <- readRDS(file.path(TXT_OBJ, "txt_metadata_raw.rds"))
mat_raw  <- readRDS(file.path(TXT_OBJ, "txt_matrix_raw.rds"))

shared_ids <- meta_raw[, .(n = uniqueN(source_lab)), by = strain_id][n > 1]$strain_id
keep <- meta_raw$strain_id %in% shared_ids
meta_sh <- meta_raw[keep]; mat_sh <- mat_raw[keep, , drop = FALSE]
meta_sh[, Site := fifelse(source_lab == "ALVAREZ", "Site A", "Site M")]
meta_sh[, Protocol := factor(treatment, levels = c("COL", "EXT"))]
cat(sprintf("Replicate-level: %d espectros (%d cepas)\n", nrow(mat_sh), uniqueN(meta_sh$strain_id)))

d_rep <- cosine_dist_mat(mat_sh)
pcoa_r <- cmdscale(d_rep, k = 2, eig = TRUE)
eigs_r <- pcoa_r$eig; eigs_r[eigs_r < 0] <- 0
pct_r  <- round(eigs_r[1:2] / sum(eigs_r) * 100, 1)
dtA <- data.table(PC1 = pcoa_r$points[, 1], PC2 = pcoa_r$points[, 2], Site = meta_sh$Site)
pA <- ggplot(dtA, aes(PC1, PC2, colour = Site)) +
  stat_ellipse(aes(group = Site), level = 0.90, linewidth = 0.4,
              linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 0.7, alpha = 0.75) +
  scale_colour_manual(values = c("Site A" = "#E64B35", "Site M" = "#4DBBD5")) +
  labs(title = "Replicate-level (individual spectra, raw / cosine)", tag = "A",
       x = paste0("PCoA1 (", pct_r[1], "%)"), y = paste0("PCoA2 (", pct_r[2], "%)")) +
  theme_paper(base_size = 8) +
  theme(plot.title = element_text(size = 8.5, face = "bold"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = ggplot2::margin(4, 6, 4, 4),
        legend.text = element_text(size = 8),
        legend.position = "bottom")

keys <- paste(meta_sh$strain_id, meta_sh$treatment, meta_sh$source_lab, sep = "__")
idx  <- split(seq_len(nrow(meta_sh)), keys)
mat_agg <- t(vapply(idx, function(i) colMeans(mat_sh[i, , drop = FALSE]), numeric(ncol(mat_sh))))
agg_meta <- rbindlist(lapply(names(idx), function(k) {
  i <- idx[[k]][1]
  data.table(sample = k, strain_id = meta_sh$strain_id[i], treatment = meta_sh$treatment[i],
             Site = meta_sh$Site[i])
}))
reps <- make_representations(mat_agg)

dist_fns <- list(
  "Raw / cosine"       = function(m) cosine_dist_mat(m),
  "TIC / Bray-Curtis"  = function(m) vegdist(m, method = "bray"),
  "Max / Euclidean"    = function(m) dist(m, method = "euclidean"),
  "CLR / Euclidean"    = function(m) dist(m, method = "euclidean"),
  "Binary / Jaccard"   = function(m) vegdist(m, method = "jaccard", binary = TRUE)
)
mats <- list("Raw / cosine" = reps$intensity_raw, "TIC / Bray-Curtis" = reps$intensity_tic,
            "Max / Euclidean" = reps$intensity_max, "CLR / Euclidean" = reps$intensity_clr,
            "Binary / Jaccard" = reps$binary)

perm_tab <- fread(file.path(TXT_TAB, "shared_strains_lab_effect_permanova.csv"))
perm_key <- c("Raw / cosine" = "intensity_raw", "TIC / Bray-Curtis" = "intensity_tic",
              "Max / Euclidean" = "intensity_max", "CLR / Euclidean" = "intensity_clr",
              "Binary / Jaccard" = "binary")
disp_rows <- rbindlist(lapply(names(mats), function(nm) {
  d <- dist_fns[[nm]](mats[[nm]])
  bd <- betadisper(d, group = agg_meta$Site)
  pv <- perm_tab[test == "permadisp_source_lab" & term == "Groups" &
                 representation == perm_key[[nm]], `Pr(>F)`]
  stopifnot(length(pv) == 1)
  data.table(representation = nm, Site = agg_meta$Site,
             dist_to_centroid = bd$distances,
             sig = ifelse(pv < 0.05, sprintf("p = %.3f", pv), sprintf("p = %.2f (n.s.)", pv)))
}))
rep_order <- c("Raw / cosine", "TIC / Bray-Curtis", "Max / Euclidean",
               "CLR / Euclidean", "Binary / Jaccard")
disp_rows[, representation := factor(representation, levels = rep_order)]
lab_dt <- unique(disp_rows[, .(representation, sig)])

pB <- ggplot(disp_rows, aes(Site, dist_to_centroid, fill = Site)) +
  geom_boxplot(width = 0.55, outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.3) +
  geom_text(data = lab_dt, aes(x = 1.5, y = Inf, label = sig), inherit.aes = FALSE,
            vjust = 1.6, size = 2.85) +
  scale_fill_manual(values = c("Site A" = "#E64B35", "Site M" = "#4DBBD5"), guide = "none") +
  facet_wrap(~ representation, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Distance to group centroid\n(PERMADISP)", tag = "B") +
  theme_paper(base_size = 8) +
  theme(plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = ggplot2::margin(4, 6, 4, 4),
        strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 8)) +
  coord_cartesian(clip = "off")

fig2 <- pA / pB + plot_layout(heights = c(1.15, 1))
save_fig(fig2, "fig2_final_hybrid", width = 17.4, height = 14.5, dpi = 600)
cat("Saved fig2_final_hybrid.png\n")
