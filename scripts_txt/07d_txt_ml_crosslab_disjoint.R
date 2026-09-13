source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(pROC)
})

N_ITER <- 100

metrics <- function(obs, prob, positive = "RES", negative = "SEN") {
  pred <- ifelse(prob > 0.5, positive, negative)
  sens <- mean(pred[obs == positive] == positive)
  spec <- mean(pred[obs == negative] == negative)
  auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, prob, quiet = TRUE,
                                                 levels = c(negative, positive), direction = "<"))),
                  error = function(e) NA_real_)
  data.table(balanced_accuracy = (sens + spec) / 2, sensitivity = sens,
             specificity = spec, auc = auc)
}

fit_predict <- function(X, y, tr, te, positive = "RES", negative = "SEN") {
  if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) return(NULL)
  ytr <- factor(y[tr])
  n <- min(table(ytr))
  fit <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500,
                      strata = ytr, sampsize = c(n, n))
  prob <- predict(fit, X[te, , drop = FALSE], type = "prob")[, positive]
  metrics(y[te], prob, positive, negative)
}

main <- function() {
  set.seed(42)
  meta <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X <- reps$intensity_raw
  stopifnot(identical(rownames(X), meta$sample))

  keep <- !is.na(meta$label_res_sen)
  meta <- meta[keep]; X <- X[keep, , drop = FALSE]
  y <- meta$label_res_sen

  alv <- unique(meta[source_lab == "ALVAREZ", strain_id])
  mal <- unique(meta[source_lab == "MALBRAN", strain_id])
  shared <- intersect(alv, mal)

  is_alv <- meta$source_lab == "ALVAREZ"
  is_mal <- meta$source_lab == "MALBRAN"

  out <- list()

  alv_strains <- unique(meta$strain_id[is_alv])
  for (i in seq_len(N_ITER)) {
    h <- sample(alv_strains, floor(0.8 * length(alv_strains)))
    r <- fit_predict(X, y, is_alv & meta$strain_id %in% h, is_alv & !(meta$strain_id %in% h))
    if (is.null(r)) next
    out[[length(out) + 1]] <- data.table(
      design = "intra_laboratorio", iter = i,
      n_train_cepas = length(h), n_test_cepas = length(alv_strains) - length(h),
      cepas_test_vistas_en_train = 0L, r)
  }

  r <- fit_predict(X, y, is_alv, is_mal)
  if (!is.null(r)) out[[length(out) + 1]] <- data.table(
    design = "naive_con_solape", iter = 1L,
    n_train_cepas = length(alv), n_test_cepas = length(mal),
    cepas_test_vistas_en_train = length(shared), r)

  te_shared <- is_mal & meta$strain_id %in% shared
  r <- fit_predict(X, y, is_alv, te_shared)
  if (!is.null(r)) out[[length(out) + 1]] <- data.table(
    design = "transferencia_instrumento", iter = 1L,
    n_train_cepas = length(alv), n_test_cepas = length(shared),
    cepas_test_vistas_en_train = length(shared), r)

  for (i in seq_len(N_ITER)) {
    half <- sample(shared, floor(length(shared) / 2))
    other <- setdiff(shared, half)
    tr <- is_alv & !(meta$strain_id %in% other)
    te <- is_mal & !(meta$strain_id %in% half)
    stopifnot(length(intersect(meta$strain_id[tr], meta$strain_id[te])) == 0)
    r <- fit_predict(X, y, tr, te)
    if (is.null(r)) next
    out[[length(out) + 1]] <- data.table(
      design = "disjunto_por_cepa", iter = i,
      n_train_cepas = uniqueN(meta$strain_id[tr]),
      n_test_cepas = uniqueN(meta$strain_id[te]),
      cepas_test_vistas_en_train = 0L, r)
  }

  res <- rbindlist(out, fill = TRUE)
  save_table(res, "crosslab_disjoint_iterations")

  summ <- res[, .(n_iter = .N,
                  cepas_train = round(mean(n_train_cepas)),
                  cepas_test = round(mean(n_test_cepas)),
                  solape_cepas = round(mean(cepas_test_vistas_en_train)),
                  BA = round(mean(balanced_accuracy, na.rm = TRUE), 3),
                  BA_sd = round(sd(balanced_accuracy, na.rm = TRUE), 3),
                  AUC = round(mean(auc, na.rm = TRUE), 3),
                  sens = round(mean(sensitivity, na.rm = TRUE), 3),
                  spec = round(mean(specificity, na.rm = TRUE), 3)), by = design]
  save_table(summ, "crosslab_disjoint_summary")
  print(summ)
  cat("\nSalidas en:", TXT_MODULAR_TAB, "\n")
}

main()
