source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(pROC)
  library(vegan)
})

N_ITER <- 100
OUT <- file.path(dirname(TXT_MODULAR_OUT), "SITE_M")
dir.create(OUT, showWarnings = FALSE)
PANEL <- c(2416, 2692, 3180, 4154)

cosine_d <- function(x) {
  x <- x / sqrt(rowSums(x^2))
  as.dist(pmax(1 - tcrossprod(x), 0))
}

prevalence_filter <- function(X, y) {
  p <- sapply(unique(y), function(k) colMeans(X[y == k, , drop = FALSE] > 0))
  X[, apply(p, 1, max) >= 0.05, drop = FALSE]
}

cliffs <- function(a, b) {
  u <- wilcox.test(a, b, exact = FALSE)$statistic
  as.numeric(2 * u / (length(a) * length(b)) - 1)
}

permanova_site <- function(X, m, label) {
  Xf <- prevalence_filter(X, m$label_res_sen)
  Xf <- Xf[rowSums(Xf) > 0, , drop = FALSE]
  mm <- as.data.frame(m[match(rownames(Xf), m$sample)])
  d <- cosine_d(Xf)
  fit <- adonis2(d ~ label_res_sen + treatment, data = mm, permutations = 999, by = "margin")
  disp <- permutest(betadisper(d, mm$label_res_sen), permutations = 999)
  data.table(site = label, n = nrow(mm), bins = ncol(Xf), term = rownames(fit)[1:2],
             R2 = round(fit$R2[1:2], 3), F = round(fit$F[1:2], 1), p = fit$`Pr(>F)`[1:2],
             permadisp_p = disp$tab$`Pr(>F)`[1])
}

classify_site <- function(X, m, label) {
  strains <- unique(m$strain_id)
  y <- m$label_res_sen
  res <- rbindlist(lapply(seq_len(N_ITER), function(i) {
    h <- sample(strains, floor(0.8 * length(strains)))
    tr <- m$strain_id %in% h; te <- !tr
    if (uniqueN(y[te]) < 2) return(NULL)
    ytr <- factor(y[tr]); n <- min(table(ytr))
    fit <- randomForest(x = X[tr, ], y = ytr, ntree = 500, strata = ytr, sampsize = c(n, n))
    prob <- predict(fit, X[te, ], type = "prob")[, "RES"]
    pred <- ifelse(prob > 0.5, "RES", "SEN")
    sens <- mean(pred[y[te] == "RES"] == "RES"); spec <- mean(pred[y[te] == "SEN"] == "SEN")
    data.table(iter = i, BA = (sens + spec) / 2, sens = sens, spec = spec,
               AUC = as.numeric(auc(roc(y[te], prob, quiet = TRUE, levels = c("SEN", "RES"), direction = "<"))))
  }))
  res[, .(site = label, n_iter = .N, BA = round(mean(BA), 3), BA_sd = round(sd(BA), 3),
          BA_lo = round(quantile(BA, .025), 3), BA_hi = round(quantile(BA, .975), 3),
          AUC = round(mean(AUC), 3), AUC_sd = round(sd(AUC), 3),
          sens = round(mean(sens), 3), spec = round(mean(spec), 3))]
}

bin_tests <- function(X, y) {
  keep <- which(colMeans(X > 0) >= 0.05)
  r <- rbindlist(lapply(keep, function(j) {
    a <- X[y == "RES", j]; b <- X[y == "SEN", j]
    data.table(bin = as.numeric(colnames(X)[j]),
               p = wilcox.test(a, b, exact = FALSE)$p.value, delta = cliffs(a, b),
               prev_res = mean(a > 0), prev_sen = mean(b > 0))
  }))
  r[, fdr := p.adjust(p, "BH")]
  r
}

main <- function() {
  set.seed(42)
  meta <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  X <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))$intensity_raw
  keep <- !is.na(meta$label_res_sen)
  meta <- meta[keep]; X <- X[keep, , drop = FALSE]
  A <- meta$source_lab == "ALVAREZ"; M <- meta$source_lab == "MALBRAN"

  comp <- meta[, .(obs = .N, cepas = uniqueN(strain_id)), by = .(source_lab, label_res_sen)]
  print(comp); fwrite(comp, file.path(OUT, "composicion.csv"))

  perm <- rbind(permanova_site(X[A, ], meta[A], "A"), permanova_site(X[M, ], meta[M], "M"))
  print(perm); fwrite(perm, file.path(OUT, "permanova_por_sitio.csv"))

  cl <- rbind(classify_site(X[A, ], meta[A], "A"), classify_site(X[M, ], meta[M], "M"))
  print(cl); fwrite(cl, file.path(OUT, "clasificacion_por_sitio.csv"))

  ba <- bin_tests(X[A, ], meta$label_res_sen[A]); bm <- bin_tests(X[M, ], meta$label_res_sen[M])
  j <- merge(ba, bm, by = "bin", suffixes = c("_A", "_M"))
  fwrite(j, file.path(OUT, "biomarcadores_por_sitio.csv"))
  sig <- j[fdr_A < 0.05]
  conc <- data.table(bins_comunes = nrow(j), sig_A = nrow(sig), sig_M = j[fdr_M < 0.05, .N],
                     sig_ambos_mismo_signo = sig[fdr_M < 0.05 & sign(delta_A) == sign(delta_M), .N],
                     sig_ambos_signo_opuesto = sig[fdr_M < 0.05 & sign(delta_A) != sign(delta_M), .N],
                     spearman_delta = round(cor(j$delta_A, j$delta_M, method = "spearman"), 3))
  print(conc); fwrite(conc, file.path(OUT, "concordancia_biomarcadores.csv"))
  pan <- j[bin %in% PANEL, .(bin, delta_A = round(delta_A, 3), fdr_A = signif(fdr_A, 2),
                             delta_M = round(delta_M, 3), fdr_M = signif(fdr_M, 2),
                             prev_res_M = round(prev_res_M, 2), prev_sen_M = round(prev_sen_M, 2))]
  print(pan); fwrite(pan, file.path(OUT, "panel_crkp_cskp_por_sitio.csv"))
}

main()
