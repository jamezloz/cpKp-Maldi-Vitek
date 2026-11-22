# Sensitivity of the main results to strain codes with suffixes (Sections 2.1, 2.9).
# Eleven strains have a second entry under a suffixed code (e.g. HP01 and HP01.1).
# Internal classification, the cross-site designs of 07d and the ST307 control of 19
# are repeated with partitions blocked by code and by base code (same seeds).
source("scripts_txt/00_shared_utils.R")
suppressPackageStartupMessages({ library(data.table); library(randomForest); library(parallel) })
RNGkind("L'Ecuyer-CMRG")
N_ITER <- 100; CORES <- 6
OUT <- file.path(dirname(TXT_MODULAR_OUT), "SUFFIX_SENSITIVITY"); dir.create(OUT, showWarnings = FALSE)
base_code <- function(s) sub("\\.[0-9.]+$", "", s)

prev_filter <- function(X, y, p = 0.05) {
  keep <- apply(X > 0, 2, function(v) max(tapply(v, y, mean)) >= p)
  X[, keep, drop = FALSE]
}

metrics <- function(yte, prob, pos, neg) {
  pred <- ifelse(prob >= 0.5, pos, neg)
  sens <- mean(pred[yte == pos] == pos); spec <- mean(pred[yte == neg] == neg)
  auc <- as.numeric(pROC::auc(pROC::roc(as.numeric(yte == pos), prob, quiet = TRUE, direction = "<")))
  data.table(BA = (sens + spec) / 2, AUC = auc, sens = sens, spec = spec)
}

cv_internal <- function(X, y, g, pos, neg, seed) {
  ug <- unique(g)
  set.seed(seed)
  res <- mclapply(seq_len(N_ITER), function(i) {
    set.seed(seed + i)
    trg <- sample(ug, floor(0.8 * length(ug))); tr <- g %in% trg; te <- !tr
    if (uniqueN(y[tr]) < 2 || uniqueN(y[te]) < 2) return(NULL)
    p <- apply(X[tr, , drop = FALSE], 2, function(v) tryCatch(wilcox.test(v ~ y[tr])$p.value, error = function(e) NA))
    q <- p.adjust(p, "BH"); sel <- which(!is.na(q) & q < 0.05)
    if (length(sel) < 2) sel <- order(p)[1:10]
    fit <- randomForest(x = X[tr, sel, drop = FALSE], y = factor(y[tr]), ntree = 500)
    metrics(y[te], predict(fit, X[te, sel, drop = FALSE], type = "prob")[, pos], pos, neg)
  }, mc.cores = CORES)
  rbindlist(res)
}

fit_bal <- function(X, y, tr, te) {
  if (uniqueN(y[tr]) < 2 || uniqueN(y[te]) < 2) return(NULL)
  ytr <- factor(y[tr]); n <- min(table(ytr))
  fit <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500, strata = ytr, sampsize = c(n, n))
  metrics(y[te], predict(fit, X[te, , drop = FALSE], type = "prob")[, "RES"], "RES", "SEN")
}

crosslab <- function(X, y, meta, g, seed) {
  is_a <- meta$source_lab == "ALVAREZ"; is_m <- meta$source_lab == "MALBRAN"
  ga <- unique(g[is_a]); shared <- intersect(ga, unique(g[is_m]))
  intra <- rbindlist(mclapply(seq_len(N_ITER), function(i) {
    set.seed(seed + i); h <- sample(ga, floor(0.8 * length(ga)))
    fit_bal(X, y, is_a & g %in% h, is_a & !(g %in% h))
  }, mc.cores = CORES))
  disj <- rbindlist(mclapply(seq_len(N_ITER), function(i) {
    set.seed(seed + 1000 + i); half <- sample(shared, floor(length(shared) / 2)); other <- setdiff(shared, half)
    tr <- is_a & !(g %in% other); te <- is_m & !(g %in% half)
    stopifnot(length(intersect(g[tr], g[te])) == 0)
    fit_bal(X, y, tr, te)
  }, mc.cores = CORES))
  list(intra = intra, disjoint = disj)
}

summ <- function(d, analisis, agrup) d[, .(analysis = analisis, grouping = agrup, n_iter = .N,
  BA = round(mean(BA), 3), BA_sd = round(sd(BA), 3), AUC = round(mean(AUC), 3),
  sens = round(mean(sens), 3), spec = round(mean(spec), 3))]

main <- function() {
  tm <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  tx <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))$intensity_raw
  stopifnot(identical(rownames(tx), tm$sample))
  rd <- readRDS(file.path(dirname(TXT_MODULAR_OBJ), "../RAW_modular/objects/raw_peak_strain_treatment_dataset.rds"))
  rm_ <- as.data.table(rd$metadata); rx <- rd$matrix
  out <- list()
  for (agr in c("code", "base_code")) {
    gf <- if (agr == "code") identity else base_code
    for (ds in c("peak-list", "RAW")) {
      M <- if (ds == "RAW") rm_ else tm; X <- if (ds == "RAW") rx else tx
      k <- !is.na(M$label_res_sen); y <- M$label_res_sen[k]
      Xk <- prev_filter(X[k, ], y)
      out[[length(out) + 1]] <- summ(cv_internal(Xk, y, gf(M$strain_id[k]), "RES", "SEN", 42), paste(ds, "CRKp/CSKp internal"), agr)
      k <- M$mechanism %in% c("KPC", "NDM"); y <- M$mechanism[k]
      Xk <- prev_filter(X[k, ], y)
      out[[length(out) + 1]] <- summ(cv_internal(Xk, y, gf(M$strain_id[k]), "KPC", "NDM", 42), paste(ds, "KPC/NDM internal"), agr)
    }
    k <- !is.na(tm$label_res_sen)
    cl <- crosslab(tx[k, ], tm$label_res_sen[k], tm[k], gf(tm$strain_id[k]), 42)
    out[[length(out) + 1]] <- summ(cl$intra, "peak-list intra-site A (07d)", agr)
    out[[length(out) + 1]] <- summ(cl$disjoint, "peak-list strain-disjoint A->M (07d)", agr)
    k <- tm$st %in% "307" & !is.na(tm$label_res_sen)
    out[[length(out) + 1]] <- summ(cv_internal(tx[k, ], tm$label_res_sen[k], gf(tm$strain_id[k]), "RES", "SEN", 42), "ST307 CRKp/CSKp (19)", agr)
    cat("done:", agr, "\n")
  }
  res <- rbindlist(out)
  fwrite(res, file.path(OUT, "suffix_sensitivity.csv"))
  print(res)
}
main()
