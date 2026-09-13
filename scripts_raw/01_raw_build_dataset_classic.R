source('scripts_raw/00_shared_utils_raw.R')

main <- function() {

  meta <- build_raw_metadata('primary')

  baseline <- 'tophat'
  specs <- lapply(meta$path, preprocess_raw_spectrum, baseline = baseline, normalize = FALSE)
  names(specs) <- meta$sample

  peak_obj <- build_peak_matrix_raw(specs)
  peak_mat <- peak_obj$matrix[meta$sample, , drop = FALSE]

  peak_det <- peak_obj$detected[meta$sample, , drop = FALSE]
  reps <- make_raw_peak_representations(peak_mat, peak_det)

  save_raw_object(meta, 'raw_metadata')
  save_raw_object(specs, 'raw_specs_preprocessed')
  save_raw_object(peak_obj$peaks, 'raw_peaks_detected')
  save_raw_object(peak_obj$peaks_binned, 'raw_peaks_binned')
  save_raw_object(peak_mat, 'raw_peak_matrix')
  save_raw_object(peak_det, 'raw_peak_detected')
  save_raw_object(reps, 'raw_peak_representations')

  save_raw_table(meta, 'raw_metadata')
  save_raw_table(as.data.table(peak_mat, keep.rownames = 'sample'), 'raw_peak_matrix')
  save_raw_table(data.table(baseline = baseline, smoothing = 'SavitzkyGolay_hws10',
                            normalization = 'ninguna (TIC pasa a ser representacion)',
                            peak_detection = 'adaptive_SNR', binning_tolerance = 0.003,
                            criterio_seleccion = 'consistencia entre replicas, ciega a la etiqueta'),
                 'raw_classic_pipeline_settings')

  agg_rep <- aggregate_raw_by_group(peak_mat, meta, 'replicate')
  agg_strain <- aggregate_raw_by_group(peak_mat, meta, 'strain_mean')
  agg_st <- aggregate_raw_by_group(peak_mat, meta, 'strain_treatment')
  agg_rep$detected <- peak_det
  agg_strain$detected <- (aggregate_raw_by_group(peak_det, meta, 'strain_mean')$matrix >= 0.5) * 1
  agg_st$detected <- (aggregate_raw_by_group(peak_det, meta, 'strain_treatment')$matrix >= 0.5) * 1
  save_raw_object(agg_rep, 'raw_peak_replicate_dataset')
  save_raw_object(agg_strain, 'raw_peak_strain_mean_dataset')
  save_raw_object(agg_st, 'raw_peak_strain_treatment_dataset')
  save_raw_table(agg_rep$metadata, 'raw_peak_replicate_metadata')
  save_raw_table(agg_strain$metadata, 'raw_peak_strain_mean_metadata')
  save_raw_table(agg_st$metadata, 'raw_peak_strain_treatment_metadata')

  save_raw_table(as.data.table(agg_st$matrix, keep.rownames = 'sample'),
                 'raw_peak_strain_treatment_matrix')

  cat('Saved RAW classic dataset objects to:', RAW_MODULAR_OUT, '\n')
  cat('RAW spectra:', nrow(meta), '\n')
  cat('RAW strain_mean samples:', nrow(agg_strain$metadata), '\n')
  cat('RAW strain_treatment samples:', nrow(agg_st$metadata), '\n')
}

main()
