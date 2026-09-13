source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
})

PERMUTATIONS <- 999

tidy_fit <- function(fit, representation, test) {
  data.table(representation = representation, test = test,
             term = rownames(fit), as.data.table(fit, keep.rownames = NULL))
}

run_models <- function(name, d, meta) {
  safe_require("vegan")
  out <- list()

  out[[length(out) + 1]] <- tidy_fit(
    vegan::adonis2(d ~ label_res_sen + mechanism + source_lab + treatment,
                   data = meta, permutations = PERMUTATIONS),
    name, "composite")

  out[[length(out) + 1]] <- tidy_fit(
    vegan::adonis2(d ~ label_res_sen + source_lab + treatment,
                   data = meta, permutations = PERMUTATIONS, by = "margin"),
    name, "margin_res_sen")

  res_only <- which(meta$label_res_sen == "RES")
  if (length(res_only) > 20 && uniqueN(meta$mechanism[res_only]) > 1) {
    d_res <- as.dist(as.matrix(d)[res_only, res_only, drop = FALSE])
    meta_res <- droplevels(as.data.frame(meta[res_only]))
    fit <- try(vegan::adonis2(d_res ~ mechanism + source_lab + treatment,
                              data = meta_res, permutations = PERMUTATIONS, by = "margin"),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
      out[[length(out) + 1]] <- tidy_fit(fit, name, "margin_mechanism")
    }
  }

  out[[length(out) + 1]] <- tidy_fit(
    vegan::adonis2(d ~ label_res_sen + mechanism + source_lab,
                   data = meta, permutations = PERMUTATIONS, strata = meta$treatment),
    name, "stratified_treatment")

  keep_kn <- meta$mechanism %in% c("KPC", "NDM")
  if (sum(keep_kn) > 20 && uniqueN(meta$mechanism[keep_kn]) > 1) {
    d_kn <- as.dist(as.matrix(d)[keep_kn, keep_kn, drop = FALSE])
    meta_kn <- droplevels(as.data.frame(meta[keep_kn, ]))
    fit_kn <- try(vegan::adonis2(d_kn ~ mechanism + source_lab + treatment,
                                 data = meta_kn, permutations = PERMUTATIONS, by = "margin"),
                 silent = TRUE)
    if (!inherits(fit_kn, "try-error")) {
      out[[length(out) + 1]] <- tidy_fit(fit_kn, name, "margin_kpc_ndm")
    }
  }

  rbindlist(out, fill = TRUE)
}

main <- function() {
  set.seed(42)
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  dist_files <- list.files(TXT_MODULAR_OBJ, pattern = "^distance_.*\\.rds$", full.names = TRUE)

  res <- rbindlist(lapply(dist_files, function(f) {
    nm <- sub("^distance_|\\.rds$", "", basename(f))
    cat("PERMANOVA:", nm, "\n")
    run_models(nm, readRDS(f), meta)
  }), fill = TRUE)

  save_table(res, "permanova_models")

  cat("\n=== R2 parcial por factor (margin_res_sen) ===\n")
  print(dcast(res[test == "margin_res_sen" & term %in% c("label_res_sen", "source_lab", "treatment")],
              representation ~ term, value.var = "R2"))
  cat("\n=== composite (modelo completo) ===\n")
  print(res[test == "composite" & term == "Model", .(representation, R2 = round(R2, 4), p = `Pr(>F)`)])
  cat("\n=== tipo de carbapenemasa, solo resistentes ===\n")
  print(res[test == "margin_mechanism" & term == "mechanism",
            .(representation, R2 = round(R2, 4), p = `Pr(>F)`)])
}

main()
