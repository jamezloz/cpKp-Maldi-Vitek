suppressPackageStartupMessages({
  library(data.table)
})

source("00_paper1_metadata.R")

TXT_MODULAR_OUT <- "resultados/TXT_modular"
TXT_MODULAR_OBJ <- file.path(TXT_MODULAR_OUT, "objects")
TXT_MODULAR_TAB <- file.path(TXT_MODULAR_OUT, "tables")
TXT_MODULAR_FIG <- file.path(TXT_MODULAR_OUT, "figures")
dir.create(TXT_MODULAR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TXT_MODULAR_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(TXT_MODULAR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(TXT_MODULAR_FIG, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

safe_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Package '%s' is required.", pkg))
}

build_txt_metadata <- function(paths, source_lab) {
  build_metadata(paths, source_lab)[!is.na(label_res_sen)]
}

read_peak_txt <- function(path) {
  dt <- tryCatch(fread(path, header = FALSE), error = function(e) NULL)
  if (is.null(dt) || ncol(dt) < 2) return(NULL)
  setnames(dt, c("mz", "intensity"))
  dt <- dt[!is.na(mz) & !is.na(intensity)]
  dt[]
}

bin_txt_spectrum <- function(dt, bin_width = 2) {
  if (is.null(dt) || nrow(dt) == 0) return(data.table(bin = numeric(), intensity = numeric()))
  tmp <- copy(dt)
  tmp[, bin := round(mz / bin_width) * bin_width]
  tmp[, .(intensity = max(intensity, na.rm = TRUE)), by = bin][order(bin)]
}

build_txt_matrix <- function(metadata, bin_width = 2) {
  spectra <- lapply(metadata$path, function(p) bin_txt_spectrum(read_peak_txt(p), bin_width = bin_width))
  all_bins <- sort(unique(unlist(lapply(spectra, function(x) x$bin))))
  mat <- matrix(0, nrow = nrow(metadata), ncol = length(all_bins))
  colnames(mat) <- as.character(all_bins)
  rownames(mat) <- metadata$sample
  for (i in seq_along(spectra)) {
    sp <- spectra[[i]]
    if (nrow(sp) == 0) next
    idx <- match(as.character(sp$bin), colnames(mat))
    mat[i, idx] <- sp$intensity
  }
  mat
}

aggregate_replicates_by_strain_treatment <- function(mat, metadata) {
  first_or_na <- function(x) {
    x <- unique(x[!is.na(x)])
    if (length(x) == 0) return(NA_character_)
    as.character(x[1])
  }
  mode_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }

  keys <- paste(metadata$strain_id, metadata$treatment, metadata$source_lab, sep = "__")
  split_idx <- split(seq_len(nrow(metadata)), keys)
  agg <- t(vapply(split_idx, function(idx) colMeans(mat[idx, , drop = FALSE]), numeric(ncol(mat))))
  meta <- rbindlist(lapply(split_idx, function(idx) {
    data.table(
      sample = paste(first_or_na(metadata$strain_id[idx]), first_or_na(metadata$treatment[idx]),
                     first_or_na(metadata$source_lab[idx]), sep = "__"),
      strain_id = first_or_na(metadata$strain_id[idx]),
      treatment = first_or_na(metadata$treatment[idx]),
      source_lab = first_or_na(metadata$source_lab[idx]),
      st = first_or_na(metadata$st[idx]),
      mechanism = mode_or_na(metadata$mechanism[idx]),
      mechanism_raw = mode_or_na(metadata$mechanism_raw[idx]),
      label_res_sen = mode_or_na(metadata$label_res_sen[idx])
    )
  }), fill = TRUE)
  rownames(agg) <- meta$sample
  list(matrix = agg, metadata = meta)
}

aggregate_replicates_by_strain <- function(mat, metadata) {
  first_or_na <- function(x) {
    x <- unique(x[!is.na(x)])
    if (length(x) == 0) return(NA_character_)
    as.character(x[1])
  }
  mode_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }
  keys <- paste(metadata$strain_id, metadata$source_lab, sep = "__")
  split_idx <- split(seq_len(nrow(metadata)), keys)
  agg <- t(vapply(split_idx, function(idx) colMeans(mat[idx, , drop = FALSE]), numeric(ncol(mat))))
  meta <- rbindlist(lapply(split_idx, function(idx) {
    data.table(
      sample = paste(first_or_na(metadata$strain_id[idx]), first_or_na(metadata$source_lab[idx]), sep = "__"),
      strain_id = first_or_na(metadata$strain_id[idx]),
      source_lab = first_or_na(metadata$source_lab[idx]),
      n_treatments = length(unique(metadata$treatment[idx])),
      treatments = paste(sort(unique(metadata$treatment[idx])), collapse = "+"),
      st = first_or_na(metadata$st[idx]),
      mechanism = mode_or_na(metadata$mechanism[idx]),
      mechanism_raw = mode_or_na(metadata$mechanism_raw[idx]),
      label_res_sen = mode_or_na(metadata$label_res_sen[idx])
    )
  }), fill = TRUE)
  rownames(agg) <- meta$sample
  list(matrix = agg, metadata = meta)
}

binary_representation <- function(x) (x > 0) * 1

max_normalize <- function(x) {
  mx <- apply(x, 1, max)
  mx[mx == 0] <- 1
  x / mx
}

tic_normalize <- function(x) {
  rs <- rowSums(x)
  rs[rs == 0] <- 1
  x / rs
}

clr_transform <- function(x, pseudocount = 1e-6) {
  x <- x + pseudocount
  lx <- log(x)
  sweep(lx, 1, rowMeans(lx), "-")
}

make_representations <- function(mat) {
  list(
    binary = binary_representation(mat),
    intensity_raw = mat,
    intensity_max = max_normalize(mat),
    intensity_tic = tic_normalize(mat),
    intensity_clr = clr_transform(tic_normalize(mat))
  )
}

save_object <- function(obj, name) {
  saveRDS(obj, file.path(TXT_MODULAR_OBJ, paste0(name, ".rds")))
}

save_table <- function(dt, name) {
  out_path <- file.path(TXT_MODULAR_TAB, paste0(name, ".csv"))
  dt <- as.data.table(dt)
  if (ncol(dt) == 0) {
    fwrite(data.table(note = "empty_result"), out_path)
    return(invisible(out_path))
  }
  fwrite(dt, out_path)
  invisible(out_path)
}

save_plot <- function(plot_obj, name, width = 8, height = 6, dpi = 300) {
  out_path <- file.path(TXT_MODULAR_FIG, paste0(name, ".png"))
  ggplot2::ggsave(filename = out_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  invisible(out_path)
}

collapse_feature_windows <- function(features, window_da = 10) {
  mz <- suppressWarnings(as.numeric(features))
  ord <- order(mz)
  mz <- mz[ord]
  feats <- features[ord]
  groups <- vector("integer", length(mz))
  current_group <- 1L
  groups[1] <- current_group
  anchor <- mz[1]
  if (length(mz) > 1) {
    for (i in 2:length(mz)) {
      if (is.na(mz[i]) || is.na(anchor) || abs(mz[i] - anchor) > window_da) {
        current_group <- current_group + 1L
        anchor <- mz[i]
      }
      groups[i] <- current_group
    }
  }
  dt <- data.table(feature = feats, mz = mz, cluster_id = groups)
  dt[, `:=`(
    cluster_center = round(mean(mz, na.rm = TRUE), 0),
    cluster_min = min(mz, na.rm = TRUE),
    cluster_max = max(mz, na.rm = TRUE),
    cluster_size = .N
  ), by = cluster_id]
  dt[]
}

collapse_ranked_biomarkers <- function(dt, feature_col = "feature", score_cols = c("fdr", "p_value"), window_da = 10) {
  stopifnot(feature_col %in% names(dt))
  cols <- intersect(score_cols, names(dt))
  if (length(cols) == 0) cols <- names(dt)[1]
  tmp <- as.data.table(copy(dt))
  tmp[, feature_chr := as.character(get(feature_col))]
  clusters <- collapse_feature_windows(tmp$feature_chr, window_da = window_da)
  tmp <- merge(tmp, clusters, by.x = "feature_chr", by.y = "feature", all.x = TRUE, sort = FALSE)
  setorderv(tmp, cols, order = rep(1L, length(cols)), na.last = TRUE)
  out <- tmp[, .SD[1], by = cluster_id]
  setorderv(out, cols, order = rep(1L, length(cols)), na.last = TRUE)
  out[]
}
