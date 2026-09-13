suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

TXT_OBJ  <- "resultados/TXT_modular/objects"
TXT_TAB  <- "resultados/TXT_modular/tables"
RAW_OBJ  <- "resultados/RAW_modular/objects"
RAW_TAB  <- "resultados/RAW_modular/tables"
FIG_OUT  <- "resultados/PAPER_FIGURES"
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

COL_RES_SEN <- c(RES = "#C0392B", SEN = "#2980B9")
COL_MECH    <- c(KPC = "#E67E22", NDM = "#27AE60", SEN = "#2980B9",
                 "KPC-NDM" = "#8E44AD")
COL_LAB     <- c(ALVAREZ = "#E74C3C", MALBRAN = "#3498DB")
COL_REP     <- c(peak_binary         = "#95A5A6",
                 peak_intensity_raw  = "#C0392B",
                 peak_intensity_tic  = "#2980B9",
                 peak_intensity_max  = "#27AE60")

theme_paper <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text             = element_text(family = "sans"),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 1, colour = "black"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "#F4F4F4", colour = "grey70"),
      strip.text       = element_text(size = base_size, face = "bold"),
      plot.tag         = element_text(size = base_size + 4, face = "bold"),
      plot.margin      = ggplot2::margin(12, 14, 12, 12),
      panel.border     = element_rect(colour = "grey50", fill = NA,
                                      linewidth = 0.6),
      legend.background = element_rect(fill = "white", colour = "grey80",
                                       linewidth = 0.3)
    )
}

save_fig <- function(p, name, width = 28, height = 20, dpi = 300) {
  path <- file.path(FIG_OUT, paste0(name, ".png"))
  ggsave(path, plot = p, width = width, height = height,
         units = "cm", dpi = dpi, bg = "white")
  cat("Saved:", path, "\n")
}

save_panel <- function(p, name, width = 15, height = 13, dpi = 300) {
  path <- file.path(FIG_OUT, paste0(name, ".png"))
  ggsave(path, plot = p, width = width, height = height,
         units = "cm", dpi = dpi, bg = "white")
  cat("  Panel saved:", path, "\n")
}
