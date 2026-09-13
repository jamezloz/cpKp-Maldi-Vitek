source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(pROC)
})

N_ITER <- 100
TRAIN_FRAC <- 0.8

load_meta_with_st <- function() {
  meta <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  if (!"st" %in% names(meta)) {

    raw <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds")))
    st_by_strain <- raw[, .(st = st[1]), by = strain_id]
    meta[, st := st_by_strain$st[match(strain_id, st_by_strain$strain_id)]]
  }
  meta[]
}

cramer_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab))
  n <- sum(tab)
  unname(sqrt(chi$statistic / (n * (min(dim(tab)) - 1))))
}

report_confounding <- function(meta) {
  s <- unique(meta[!is.na(st) & !is.na(label_res_sen), .(strain_id, st, mechanism, label_res_sen)])
  comp <- dcast(s, st ~ label_res_sen, fun.aggregate = length, value.var = "strain_id")
  comp[, total := RES + SEN]
  setorder(comp, -total)

  tab <- table(s$st, s$label_res_sen)
  chi <- suppressWarnings(chisq.test(tab))
  stats <- data.table(
    n_cepas_con_ST = nrow(s),
    n_ST = uniqueN(s$st),
    cramer_v = round(cramer_v(tab), 3),
    chisq_p = chi$p.value
  )
  save_table(comp, "clonality_st_class_composition")
  save_table(stats, "clonality_association_stats")
  cat("Composicion ST x clase:\n"); print(comp)
  cat("\nCramer V =", stats$cramer_v, " p =", format.pval(stats$chisq_p, 3), "\n")
  invisible(comp)
}

permanova_with_st <- function(meta) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NULL)

  canonical_reps <- c("intensity_raw_cosine", "intensity_tic_bray",
                      "intensity_max_euclidean", "intensity_clr_euclidean",
                      "binary_jaccard")
  dist_files <- file.path(TXT_MODULAR_OBJ, paste0("distance_", canonical_reps, ".rds"))
  stopifnot(all(file.exists(dist_files)))
  out <- list()
  for (f in dist_files) {
    nm <- sub("^distance_|\\.rds$", "", basename(f))
    d <- readRDS(f)
    keep <- !is.na(meta$st) & !is.na(meta$label_res_sen)

    st_ok <- names(which(table(meta$st[keep]) >= 10))
    keep <- keep & meta$st %in% st_ok
    if (sum(keep) < 30) next
    dm <- as.dist(as.matrix(d)[keep, keep, drop = FALSE])
    md <- droplevels(as.data.frame(meta[keep]))
    md$st <- factor(md$st)
    fit <- try(vegan::adonis2(dm ~ st + label_res_sen + treatment,
                              data = md, permutations = 999, by = "margin"), silent = TRUE)
    if (inherits(fit, "try-error")) next
    out[[length(out) + 1]] <- data.table(representation = nm, term = rownames(fit),
                                         as.data.table(fit, keep.rownames = NULL))
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out, fill = TRUE)
  save_table(res, "clonality_permanova_st_covariate")
  cat("\nPERMANOVA con ST como covariable (by=margin):\n")
  print(res[term %in% c("st", "label_res_sen", "treatment"),
            .(representation, term, R2 = round(R2, 4), p = `Pr(>F)`)])
  invisible(res)
}

grouped_split <- function(groups, train_frac = TRAIN_FRAC) {
  ug <- unique(groups)
  tr_g <- sample(ug, size = max(2, floor(train_frac * length(ug))))
  groups %in% tr_g
}

run_rf_iterations <- function(X, y, groups, n_iter = N_ITER, label = "") {
  res <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    tr <- grouped_split(groups)
    if (length(unique(y[tr])) < 2 || length(unique(y[!tr])) < 2) next
    fit <- randomForest(x = X[tr, , drop = FALSE], y = factor(y[tr]), ntree = 500)
    prob <- predict(fit, X[!tr, , drop = FALSE], type = "prob")[, "RES"]
    pred <- ifelse(prob > 0.5, "RES", "SEN")
    obs <- y[!tr]
    sens <- mean(pred[obs == "RES"] == "RES")
    spec <- mean(pred[obs == "SEN"] == "SEN")
    auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, prob, quiet = TRUE, levels = c("SEN", "RES"), direction = "<"))),
                    error = function(e) NA_real_)
    res[[i]] <- data.table(iter = i, analysis = label,
                           balanced_accuracy = (sens + spec) / 2,
                           sensitivity = sens, specificity = spec, auc = auc)
  }
  rbindlist(res, fill = TRUE)
}

within_clone_classifier <- function(meta) {
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X_all <- reps$intensity_raw
  stopifnot(identical(rownames(X_all), meta$sample))

  out <- list()

  keep <- !is.na(meta$label_res_sen)
  out$all <- run_rf_iterations(X_all[keep, , drop = FALSE], meta$label_res_sen[keep],
                               meta$strain_id[keep], label = "todos_los_ST")

  k307 <- keep & !is.na(meta$st) & meta$st == "307"
  cat("\nST307: muestras =", sum(k307),
      " cepas =", uniqueN(meta$strain_id[k307]),
      " RES/SEN =", paste(table(meta$label_res_sen[k307]), collapse = "/"), "\n")
  if (uniqueN(meta$strain_id[k307]) >= 20) {
    out$st307 <- run_rf_iterations(X_all[k307, , drop = FALSE], meta$label_res_sen[k307],
                                   meta$strain_id[k307], label = "solo_ST307")
  }

  k307_alv <- k307 & meta$source_lab == "ALVAREZ"
  cat("ST307 solo ALVAREZ: muestras =", sum(k307_alv),
      " cepas =", uniqueN(meta$strain_id[k307_alv]),
      " RES/SEN =", paste(table(meta$label_res_sen[k307_alv]), collapse = "/"), "\n")
  if (uniqueN(meta$strain_id[k307_alv]) >= 20) {
    out$st307_alv <- run_rf_iterations(X_all[k307_alv, , drop = FALSE],
                                       meta$label_res_sen[k307_alv],
                                       meta$strain_id[k307_alv], label = "ST307_solo_Alvarez")
  }

  res <- rbindlist(out, fill = TRUE)
  save_table(res, "clonality_within_clone_iterations")

  summ <- res[, .(n_iter = .N,
                  BA = round(mean(balanced_accuracy, na.rm = TRUE), 3),
                  BA_sd = round(sd(balanced_accuracy, na.rm = TRUE), 3),
                  AUC = round(mean(auc, na.rm = TRUE), 3),
                  AUC_sd = round(sd(auc, na.rm = TRUE), 3),
                  sens = round(mean(sensitivity, na.rm = TRUE), 3),
                  spec = round(mean(specificity, na.rm = TRUE), 3)),
              by = analysis]
  save_table(summ, "clonality_within_clone_summary")
  cat("\nClasificador RES/SEN dentro del mismo clon:\n"); print(summ)
  invisible(summ)
}

run_rf_mech <- function(X, y, groups, n_iter = N_ITER, label = "") {
  res <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    tr <- grouped_split(groups)
    if (length(unique(y[tr])) < 2 || length(unique(y[!tr])) < 2) next
    fit <- randomForest(x = X[tr, , drop = FALSE], y = factor(y[tr]), ntree = 500)
    prob <- predict(fit, X[!tr, , drop = FALSE], type = "prob")[, "KPC"]
    pred <- ifelse(prob > 0.5, "KPC", "NDM"); obs <- y[!tr]
    sens <- mean(pred[obs == "KPC"] == "KPC"); spec <- mean(pred[obs == "NDM"] == "NDM")
    auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, prob, quiet = TRUE, levels = c("NDM", "KPC"), direction = "<"))),
                    error = function(e) NA_real_)
    res[[i]] <- data.table(iter = i, analysis = label,
                           balanced_accuracy = (sens + spec) / 2,
                           sensitivity = sens, specificity = spec, auc = auc)
  }
  rbindlist(res, fill = TRUE)
}

within_clone_mechanism <- function(meta) {
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X_all <- reps$intensity_raw
  base <- !is.na(meta$mechanism) & meta$mechanism %in% c("KPC", "NDM")
  out <- list()
  out$all <- run_rf_mech(X_all[base, , drop = FALSE], meta$mechanism[base],
                         meta$strain_id[base], label = "todos_los_ST")
  for (stv in c("258", "11", "307")) {
    k <- base & !is.na(meta$st) & meta$st == stv
    cat(sprintf("\nST%s KPC/NDM: muestras = %d cepas = %d KPC/NDM = %s\n", stv,
                sum(k), uniqueN(meta$strain_id[k]),
                paste(table(meta$mechanism[k]), collapse = "/")))
    if (uniqueN(meta$strain_id[k]) >= 15)
      out[[paste0("ST", stv)]] <- run_rf_mech(X_all[k, , drop = FALSE], meta$mechanism[k],
                                              meta$strain_id[k], label = paste0("solo_ST", stv))
  }
  res <- rbindlist(out, fill = TRUE)
  save_table(res, "clonality_mechanism_within_clone_iterations")
  summ <- res[, .(n_iter = .N, BA = round(mean(balanced_accuracy, na.rm = TRUE), 3),
                  BA_sd = round(sd(balanced_accuracy, na.rm = TRUE), 3),
                  AUC = round(mean(auc, na.rm = TRUE), 3),
                  sens = round(mean(sensitivity, na.rm = TRUE), 3),
                  spec = round(mean(specificity, na.rm = TRUE), 3)), by = analysis]
  save_table(summ, "clonality_mechanism_within_clone_summary")
  cat("\nClasificador KPC/NDM dentro del mismo clon:\n"); print(summ)
}

main <- function() {
  set.seed(42)
  meta <- load_meta_with_st()
  report_confounding(meta)
  permanova_with_st(meta)
  within_clone_classifier(meta)
  within_clone_mechanism(meta)
  cat("\nSalidas en:", TXT_MODULAR_TAB, "\n")
}

main()
