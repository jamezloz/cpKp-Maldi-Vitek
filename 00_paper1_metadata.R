suppressPackageStartupMessages({
  library(data.table)
})

DATA_BASE <- "PATH_TO_SPECTRA"

TXT_ALVAREZ_DIR <- file.path(DATA_BASE, "peaklists_siteA")
TXT_MALBRAN_DIR <- file.path(DATA_BASE, "peaklists_siteM")
RAW_ALVAREZ_DIR <- file.path(DATA_BASE, "mzml_siteA")
RAW_MALBRAN_DIR <- file.path(DATA_BASE, "mzml_siteM")


parse_spectrum_names <- function(paths) {
  stem <- sub("\\.[^.]+$", "", basename(paths))

  parts <- tstrsplit(stem, "_", fixed = TRUE)
  if (length(parts) != 4L) {
    bad <- basename(paths)[lengths(strsplit(stem, "_", fixed = TRUE)) != 4L]
    stop(sprintf("Nombres fuera del formato canonico (%d): %s",
                 length(bad), paste(utils::head(bad, 5), collapse = ", ")))
  }

  st_rep <- parts[[4]]
  st     <- toupper(sub("^st", "", sub("-[0-9]+$", "", st_rep)))
  rep    <- sub("^.*-", "", st_rep)
  st[st == "NA"] <- NA_character_

  data.table(
    path          = paths,
    sample        = basename(paths),
    strain_id     = toupper(parts[[1]]),
    treatment     = toupper(parts[[2]]),
    mechanism_raw = toupper(parts[[3]]),
    st            = st,
    replicate     = rep
  )
}

collapse_mechanism <- function(mech) {
  fifelse(mech == "KPC-OXA48", "KPC",
    fifelse(mech == "NDM-OXA48", "NDM",
      fifelse(mech == "NDM-KPC", "KPC-NDM", mech)))
}

label_res_sen <- function(mech) {
  fifelse(is.na(mech) | mech == "OXA48", NA_character_,
          fifelse(mech == "SEN", "SEN", "RES"))
}

build_metadata <- function(paths, source_lab) {
  if (length(paths) == 0) return(NULL)
  dt <- parse_spectrum_names(paths)
  dt[, source_lab := source_lab]
  dt[, mechanism := collapse_mechanism(mechanism_raw)]
  dt[, label_res_sen := label_res_sen(mechanism)]
  dt[]
}

cosine_dist <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- 1
  sim <- (x / norms) %*% t(x / norms)
  as.dist(1 - sim)
}

jaccard_dist <- function(x) {
  x <- (as.matrix(x) > 0) * 1
  inter <- tcrossprod(x)
  sizes <- rowSums(x)
  union <- outer(sizes, sizes, "+") - inter
  as.dist(ifelse(union == 0, 0, 1 - inter / union))
}

list_spectra <- function(dir, ext) {
  list.files(dir, pattern = paste0("\\.", ext, "$"),
             full.names = TRUE, ignore.case = TRUE)
}

load_txt_metadata <- function(labs = c("ALVAREZ", "MALBRAN")) {
  out <- list()
  if ("ALVAREZ" %in% labs)
    out$alv <- build_metadata(list_spectra(TXT_ALVAREZ_DIR, "txt"), "ALVAREZ")
  if ("MALBRAN" %in% labs)
    out$mal <- build_metadata(list_spectra(TXT_MALBRAN_DIR, "txt"), "MALBRAN")
  rbindlist(out, fill = TRUE)
}

load_raw_metadata <- function(labs = c("ALVAREZ", "MALBRAN")) {
  out <- list()
  if ("ALVAREZ" %in% labs)
    out$alv <- build_metadata(list_spectra(RAW_ALVAREZ_DIR, "mzml"), "ALVAREZ")
  if ("MALBRAN" %in% labs)
    out$mal <- build_metadata(list_spectra(RAW_MALBRAN_DIR, "mzml"), "MALBRAN")
  rbindlist(out, fill = TRUE)
}

subset_dataset <- function(meta, design = c("primary", "crosslab", "pooled")) {
  design <- match.arg(design)
  if (design == "primary") return(meta[source_lab == "ALVAREZ"])
  if (design == "pooled")  return(meta)
  shared <- meta[, .(n_lab = uniqueN(source_lab)), by = strain_id][n_lab == 2, strain_id]
  meta[strain_id %in% shared]
}
