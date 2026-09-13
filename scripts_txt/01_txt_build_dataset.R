source("scripts_txt/00_shared_utils.R")

main <- function() {

  meta_raw <- rbindlist(list(
    build_txt_metadata(list_spectra(TXT_ALVAREZ_DIR, "txt"), "ALVAREZ"),
    build_txt_metadata(list_spectra(TXT_MALBRAN_DIR, "txt"), "MALBRAN")
  ), fill = TRUE)

  mat_raw <- build_txt_matrix(meta_raw, bin_width = 2)
  agg <- aggregate_replicates_by_strain_treatment(mat_raw, meta_raw)
  reps <- make_representations(agg$matrix)

  save_table(meta_raw, "txt_metadata_raw")
  save_table(as.data.table(mat_raw, keep.rownames = "sample"), "txt_matrix_raw")
  save_table(agg$metadata, "txt_metadata_aggregated")
  save_table(as.data.table(agg$matrix, keep.rownames = "sample"), "txt_matrix_aggregated")

  save_object(meta_raw, "txt_metadata_raw")
  save_object(mat_raw, "txt_matrix_raw")
  save_object(agg$metadata, "txt_metadata_aggregated")
  save_object(agg$matrix, "txt_matrix_aggregated")
  save_object(reps, "txt_representations")

  cat("Saved TXT dataset objects to:", TXT_MODULAR_OUT, "\n")
  cat("Raw spectra:", nrow(meta_raw), "\n")
  cat("Aggregated strain+treatment samples:", nrow(agg$metadata), "\n")
}

main()
