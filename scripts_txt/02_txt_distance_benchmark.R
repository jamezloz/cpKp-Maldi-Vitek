source("scripts_txt/00_shared_utils.R")

main <- function() {
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))

  dists <- list(
    binary_jaccard = jaccard_dist(reps$binary),
    intensity_raw_cosine = cosine_dist(reps$intensity_raw),
    intensity_max_cosine = cosine_dist(reps$intensity_max),
    intensity_max_euclidean = dist(reps$intensity_max),
    intensity_tic_bray = if (requireNamespace("vegan", quietly = TRUE)) vegan::vegdist(reps$intensity_tic, method = "bray") else NULL,
    intensity_clr_euclidean = dist(reps$intensity_clr)
  )

  dists <- dists[!vapply(dists, is.null, logical(1))]
  for (nm in names(dists)) save_object(dists[[nm]], paste0("distance_", nm))

  summary_dt <- data.table(
    representation = names(dists),
    metric = c("jaccard", "cosine", "cosine", "euclidean", "bray", "euclidean")[seq_along(dists)]
  )
  save_table(summary_dt, "distance_benchmark_summary")
  cat("Saved distance objects to:", TXT_MODULAR_OUT, "\n")
}

main()
