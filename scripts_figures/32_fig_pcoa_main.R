library(ape)
library(ggplot2)
library(dplyr)
library(patchwork)
library(data.table)

source("scripts_figures/00_fig_utils.R")

OBJ_TXT <- "resultados/TXT_modular/objects"
OBJ_RAW <- "resultados/RAW_modular/objects"
TAB_TXT <- "resultados/TXT_modular/tables"
TAB_RAW <- "resultados/RAW_modular/tables"

r2_txt <- function(path, subset_fun) {
  d <- data.table::fread(file.path(path))
  v <- subset_fun(d)
  stopifnot(nrow(v) >= 1)
  as.numeric(v$R2[1])
}

R2_A <- r2_txt(file.path(TAB_TXT, "permanova_models.csv"), function(d)
  d[representation == "intensity_raw_cosine.rds" & test == "margin_res_sen" &
      term == "label_res_sen"])
R2_C <- r2_txt(file.path(TAB_TXT, "permanova_summary_raw_level.csv"), function(d)
  d[representation == "intensity_raw" & test == "margin_res_sen" & term == "label_res_sen"])
R2_B <- r2_txt(file.path(TAB_RAW, "raw_classic_permanova.csv"), function(d)
  d[group_mode == "strain_treatment" & representation == "peak_intensity_raw_cosine" &
      task == "res_sen" & term == "label_res_sen"])
R2_D <- r2_txt(file.path(TAB_RAW, "raw_classic_permanova.csv"), function(d)
  d[group_mode == "replicate" & representation == "peak_intensity_raw_cosine" &
      task == "res_sen" & term == "label_res_sen"])

annot_r2 <- function(r2) sprintf("PERMANOVA resistance: R² = %.3f, p = 0.001", r2)
cat(sprintf("R2 resistencia -> A=%.3f B=%.3f C=%.3f D=%.3f\n", R2_A, R2_B, R2_C, R2_D))
OUT_DIR <- FIG_OUT
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MECH_LABELS <- c(KPC = "KPC", NDM = "NDM", "KPC-NDM" = "KPC+NDM", SEN = "CSKp")
COL_MECH    <- c(KPC = "#E69F00", NDM = "#009E73", "KPC+NDM" = "#CC79A7", CSKp = "#0072B2")
SHAPE_TREAT <- c(COL = 16, EXT = 17)

cosine_dist <- function(X) {
  norms <- sqrt(rowSums(X^2)); norms[norms == 0] <- 1e-10
  X_n <- X / norms; as.dist(1 - X_n %*% t(X_n))
}

sc_col   <- scale_colour_manual(values = COL_MECH, name = "Group",
                                limits = names(COL_MECH), drop = FALSE)
sc_fill  <- scale_fill_manual(  values = COL_MECH, limits = names(COL_MECH),
                                drop = FALSE, guide = "none")
sc_shape <- scale_shape_manual( values = SHAPE_TREAT, name = "Preparation", drop = FALSE)

build_panel <- function(df, pve, tag, title, annot,
                        xlim = NULL, ylim = NULL, pt_size = 0.9, pt_alpha = 0.80) {
  p <- ggplot(df, aes(PC1, PC2, colour = Group, shape = Treatment)) +
    stat_ellipse(aes(group = Group, fill = Group),
                 geom = "polygon", level = 0.95, alpha = 0.08,
                 linewidth = 0.3, linetype = "dashed", show.legend = FALSE) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    sc_col + sc_fill + sc_shape +
    guides(colour = guide_legend(override.aes = list(shape = 16, alpha = 1)),
           shape  = guide_legend(override.aes = list(colour = "grey25", alpha = 1))) +
    labs(title = title, tag = tag,
         x = paste0("PCoA1 (", pve[1], "%)"),
         y = paste0("PCoA2 (", pve[2], "%)")) +
    theme_paper(base_size = 8) +
    theme(legend.position       = "none",
          plot.margin           = ggplot2::margin(4, 6, 4, 4),
          plot.title            = element_text(size = 8.5, face = "bold"),
          plot.tag              = element_text(size = 11, face = "bold"),
          plot.tag.position     = c(0.02, 0.98))

  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  x_ann <- if (is.null(xlim)) Inf  else xlim[2]
  y_ann <- if (is.null(ylim)) -Inf else ylim[1]
  p + annotate("text", x = x_ann, y = y_ann, label = annot,
               hjust = 1.02, vjust = -0.4, size = 2.85, colour = "grey30")
}

pcoa_s <- readRDS(file.path(OBJ_TXT, "intensity_raw_cosine.rds_pcoa.rds"))
meta_s <- readRDS(file.path(OBJ_TXT, "txt_metadata_aggregated.rds"))
pve_sA <- round(pcoa_s$values$Eigenvalues /
                  sum(pcoa_s$values$Eigenvalues[pcoa_s$values$Eigenvalues > 0]) * 100, 1)

dfA <- as.data.frame(pcoa_s$vectors[, 1:2]) |> setNames(c("PC1","PC2")) |>
  mutate(sample = rownames(pcoa_s$vectors)) |>
  left_join(meta_s, by = "sample") |>
  mutate(Group     = factor(recode(mechanism, !!!MECH_LABELS), levels = names(COL_MECH)),
         Treatment = factor(treatment, levels = c("COL","EXT")))

pA <- build_panel(dfA, pve_sA, "A", "Peak-list — strain×protocol level",
                  annot_r2(R2_A))

dst_sB  <- readRDS(file.path(OBJ_RAW, "raw_peak_strain_treatment_dataset.rds"))
d_sB    <- cosine_dist(dst_sB$matrix)
meta_sB <- dst_sB$metadata
pcB     <- ape::pcoa(d_sB)
pve_sB  <- round(pcB$values$Eigenvalues / sum(pcB$values$Eigenvalues[pcB$values$Eigenvalues > 0]) * 100, 1)

dfB <- as.data.frame(pcB$vectors[, 1:2]) |> setNames(c("PC1","PC2")) |>
  mutate(sample = rownames(pcB$vectors)) |>
  left_join(meta_sB, by = "sample") |>
  mutate(Group     = factor(recode(mechanism, !!!MECH_LABELS), levels = names(COL_MECH)),
         Treatment = factor(treatment, levels = c("COL","EXT")))

pB <- build_panel(dfB, pve_sB, "B", "RAW dataset — strain×protocol level",
                  annot_r2(R2_B))

cat("Computing peak-list replicate cosine distance...\n")
mat_rC  <- readRDS(file.path(OBJ_TXT, "txt_matrix_raw.rds"))
meta_rC <- readRDS(file.path(OBJ_TXT, "txt_metadata_raw.rds"))
d_rC    <- cosine_dist(mat_rC)
pcC     <- ape::pcoa(d_rC)
pve_rC  <- round(pcC$values$Eigenvalues / sum(pcC$values$Eigenvalues[pcC$values$Eigenvalues > 0]) * 100, 1)

dfC <- as.data.frame(pcC$vectors[, 1:2]) |> setNames(c("PC1","PC2")) |>
  bind_cols(meta_rC[match(rownames(mat_rC), meta_rC$sample), ]) |>
  mutate(Group     = factor(recode(mechanism, !!!MECH_LABELS), levels = names(COL_MECH)),
         Treatment = factor(treatment, levels = c("COL","EXT")))

pC <- build_panel(dfC, pve_rC, "C", "Peak-list — replicate level",
                  annot_r2(R2_C),
                  pt_size = 0.5, pt_alpha = 0.55)

d_rD    <- readRDS(file.path(OBJ_RAW, "raw_classic_distances_replicate.rds"))$peak_intensity_raw_cosine
meta_rD <- readRDS(file.path(OBJ_RAW, "raw_peak_replicate_dataset.rds"))$metadata
pcD     <- ape::pcoa(d_rD)
pve_rD  <- round(pcD$values$Eigenvalues / sum(pcD$values$Eigenvalues[pcD$values$Eigenvalues > 0]) * 100, 1)

meta_rD_ord <- meta_rD[match(attr(d_rD, "Labels"), meta_rD$sample), ]
dfD <- as.data.frame(pcD$vectors[, 1:2]) |> setNames(c("PC1","PC2")) |>
  bind_cols(meta_rD_ord) |>
  mutate(Group     = factor(recode(mechanism, !!!MECH_LABELS), levels = names(COL_MECH)),
         Treatment = factor(treatment, levels = c("COL","EXT"))) |>
  filter(!is.na(Group))

pD <- build_panel(dfD, pve_rD, "D", "RAW dataset — replicate level",
                  annot_r2(R2_D),
                  pt_size = 0.5, pt_alpha = 0.55)

leg_src <- pA +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(0.35, "cm"),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8.5, face = "bold"))
legend <- cowplot::get_legend(leg_src)

grid <- wrap_plots(pA, pB, pC, pD, ncol = 2)
fig  <- cowplot::plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.08))

out_path <- file.path(OUT_DIR, "fig_pcoa_main.png")
ggsave(out_path, fig, width = 17.4, height = 15.6, units = "cm", dpi = 600, bg = "white")
cat("Saved:", out_path, "\n")
