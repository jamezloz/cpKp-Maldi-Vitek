suppressPackageStartupMessages({
  library(data.table)
})

source("00_paper1_metadata.R")
RAW_MODULAR_OUT <- "resultados/RAW_modular"
RAW_MODULAR_OBJ <- file.path(RAW_MODULAR_OUT, 'objects')
RAW_MODULAR_TAB <- file.path(RAW_MODULAR_OUT, 'tables')
RAW_MODULAR_FIG <- file.path(RAW_MODULAR_OUT, 'figures')
dir.create(RAW_MODULAR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_MODULAR_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_MODULAR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_MODULAR_FIG, recursive = TRUE, showWarnings = FALSE)

build_raw_metadata <- function(design = c("primary", "crosslab", "pooled")) {
  design <- match.arg(design)
  labs <- if (design == "primary") "ALVAREZ" else c("ALVAREZ", "MALBRAN")
  meta <- subset_dataset(load_raw_metadata(labs), design)
  meta[!is.na(label_res_sen)]
}

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Package '%s' is required.", pkg))
}

save_raw_object <- function(obj, name) {
  saveRDS(obj, file.path(RAW_MODULAR_OBJ, paste0(name, '.rds')))
}

save_raw_table <- function(dt, name) {
  out_path <- file.path(RAW_MODULAR_TAB, paste0(name, '.csv'))
  dt <- as.data.table(dt)
  if (ncol(dt) == 0) {
    fwrite(data.table(note = 'empty_result'), out_path)
    return(invisible(out_path))
  }
  fwrite(dt, out_path)
  invisible(out_path)
}

save_raw_plot <- function(plot_obj, name, width = 8, height = 6, dpi = 300) {
  out_path <- file.path(RAW_MODULAR_FIG, paste0(name, '.png'))
  ggplot2::ggsave(filename = out_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  invisible(out_path)
}

robust_q90_mad <- function(intensity) {
  q90 <- as.numeric(stats::quantile(intensity, 0.9, na.rm = TRUE))
  mad_val <- as.numeric(stats::mad(intensity, center = median(intensity, na.rm = TRUE), constant = 1, na.rm = TRUE))
  if (is.na(mad_val) || mad_val == 0) return(0)
  q90 / mad_val
}

preprocess_raw_spectrum <- function(path,
                                    baseline = c('tophat', 'snip50', 'snip25', 'snip100'),
                                    normalize = FALSE) {
  baseline <- match.arg(baseline)
  need_pkg('MALDIquant')
  need_pkg('MALDIquantForeign')
  s <- MALDIquantForeign::importMzMl(path, verbose = FALSE)[[1]]
  s <- MALDIquant::smoothIntensity(s, method = 'SavitzkyGolay', halfWindowSize = 10)
  s <- switch(baseline,
    tophat  = MALDIquant::removeBaseline(s, method = 'TopHat'),
    snip25  = MALDIquant::removeBaseline(s, method = 'SNIP', iterations = 25),
    snip50  = MALDIquant::removeBaseline(s, method = 'SNIP', iterations = 50),
    snip100 = MALDIquant::removeBaseline(s, method = 'SNIP', iterations = 100))
  if (normalize) s <- MALDIquant::calibrateIntensity(s, method = 'TIC')
  s
}

estimate_snr <- function(spec) {
  rough <- MALDIquant::detectPeaks(spec, SNR = 2)
  n <- length(rough)
  snr <- 2.5 + 0.02 * n
  snr <- max(2.5, min(8, snr))
  snr
}

build_peak_matrix_raw <- function(specs_norm) {
  peaks <- lapply(specs_norm, function(s) {
    snr <- estimate_snr(s)
    MALDIquant::detectPeaks(s, SNR = snr, halfWindowSize = 10)
  })
  peaks_binned <- MALDIquant::binPeaks(peaks, tolerance = 0.003)

  mat <- MALDIquant::intensityMatrix(peaks_binned, specs_norm)
  rownames(mat) <- names(specs_norm)

  mat_det <- MALDIquant::intensityMatrix(peaks_binned)
  rownames(mat_det) <- names(specs_norm)
  detected <- !is.na(mat_det)
  storage.mode(detected) <- "integer"

  list(matrix = mat, detected = detected, peaks = peaks, peaks_binned = peaks_binned)
}

aggregate_raw_by_group <- function(X, metadata,
                                   group_mode = c('replicate', 'strain_mean', 'strain_treatment')) {
  group_mode <- match.arg(group_mode)
  if (group_mode == 'replicate') {
    return(list(matrix = X, metadata = copy(metadata)))
  }

  if (group_mode == 'strain_mean') {
    keys <- metadata$strain_id
    sample_of <- function(g) g$strain_id[1]
  } else {
    keys <- paste(metadata$strain_id, metadata$treatment, metadata$source_lab, sep = '__')
    sample_of <- function(g) paste(g$strain_id[1], g$treatment[1], g$source_lab[1], sep = '__')
  }

  split_idx <- split(seq_len(nrow(metadata)), keys)
  Xagg <- do.call(rbind, lapply(split_idx, function(idx) colMeans(X[idx, , drop = FALSE])))
  meta <- rbindlist(lapply(split_idx, function(idx) {
    g <- metadata[idx]
    data.table(
      sample = sample_of(g),
      strain_id = g$strain_id[1],
      st = g$st[1],
      source_lab = names(sort(table(g$source_lab), decreasing = TRUE))[1],
      treatment = names(sort(table(g$treatment), decreasing = TRUE))[1],
      mechanism = names(sort(table(g$mechanism), decreasing = TRUE))[1],
      mechanism_raw = names(sort(table(g$mechanism_raw), decreasing = TRUE))[1],
      label_res_sen = names(sort(table(g$label_res_sen), decreasing = TRUE))[1]
    )
  }))
  rownames(Xagg) <- meta$sample
  list(matrix = Xagg, metadata = meta)
}

tic_normalize <- function(x) {
  rs <- rowSums(x); rs[rs == 0] <- 1; x / rs
}

max_normalize <- function(x) {
  mx <- apply(x, 1, max); mx[mx == 0] <- 1; x / mx
}

clr_transform_raw <- function(x, pseudocount = 1e-6) {
  lx <- log(x + pseudocount)
  sweep(lx, 1, rowMeans(lx), "-")
}

make_raw_peak_representations <- function(mat, detected = NULL) {
  if (is.null(detected)) {
    stop("make_raw_peak_representations necesita la matriz de deteccion: la binaria ",
         "no puede derivarse de la de intensidad (ver build_peak_matrix_raw).")
  }
  list(
    peak_binary        = detected,
    peak_intensity_raw = mat,
    peak_intensity_tic = tic_normalize(mat),
    peak_intensity_max = max_normalize(mat),
    peak_intensity_clr = clr_transform_raw(tic_normalize(mat))
  )
}

macro_balanced_accuracy <- function(truth, pred) {
  lv <- levels(truth)
  vals <- sapply(lv, function(cl) {
    idx <- truth == cl
    if (!any(idx)) return(NA_real_)
    mean(pred[idx] == cl)
  })
  mean(vals, na.rm = TRUE)
}

macro_f1 <- function(truth, pred) {
  lv <- levels(truth)
  vals <- sapply(lv, function(cl) {
    tp <- sum(truth == cl & pred == cl)
    fp <- sum(truth != cl & pred == cl)
    fn <- sum(truth == cl & pred != cl)

    if ((tp + fn) == 0 && (tp + fp) == 0) return(NA_real_)
    prec <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
  })
  mean(vals, na.rm = TRUE)
}

filter_r_like <- function(mat, y) {
  freq <- colMeans(mat > 0)
  keep1 <- freq > 0.01 & freq < 0.5
  mat1 <- mat[, keep1, drop = FALSE]
  if (ncol(mat1) < 5) return(NULL)
  vars <- apply(mat1, 2, var)
  keep2 <- vars > stats::quantile(vars, 0.75, na.rm = TRUE)
  mat2 <- mat1[, keep2, drop = FALSE]
  if (ncol(mat2) < 5) return(NULL)
  pvals <- apply(mat2, 2, function(v) stats::wilcox.test(v ~ y)$p.value)
  keep3 <- p.adjust(pvals, method = 'fdr') < 0.05
  mat3 <- mat2[, keep3, drop = FALSE]
  if (ncol(mat3) < 5) return(NULL)
  mat3
}
