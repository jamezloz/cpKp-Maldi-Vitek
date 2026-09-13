source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
})

N_BOOT <- 100
set.seed(42)

meta_raw <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds")))
mat_raw  <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))

keys <- paste(meta_raw$strain_id, meta_raw$treatment, meta_raw$source_lab, sep = "__")
idx  <- split(seq_len(nrow(meta_raw)), keys)
X <- do.call(rbind, lapply(idx, function(i) colMeans(mat_raw[i, , drop = FALSE])))
M <- rbindlist(lapply(names(idx), function(k) {
  p <- strsplit(k, "__")[[1]]
  meta_raw[idx[[k]][1], .(sample = k, strain_id, treatment, source_lab, mechanism)]
}))
stopifnot(nrow(X) == nrow(M))

ok <- M$mechanism %in% c("KPC", "NDM")
X <- X[ok, , drop = FALSE]; M <- M[ok]

shared_ids <- M[, .(n = uniqueN(source_lab)), by = strain_id][n > 1]$strain_id
cat(sprintf("Cepas KPC/NDM compartidas: %d\n", length(shared_ids)))
print(unique(M[strain_id %in% shared_ids, .(strain_id, mechanism)])[, .N, by = mechanism])

M[, sp := paste(strain_id, treatment, sep = "__")]
pairs_sp <- M[strain_id %in% shared_ids,
              .(n_sites = uniqueN(source_lab)), by = sp][n_sites == 2]$sp
cat(sprintf("Pares (cepa, protocolo) medidos en ambos sitios: %d\n", length(pairs_sp)))

train_rows <- which(!(M$strain_id %in% shared_ids) & M$source_lab == "ALVAREZ")
cat(sprintf("Entrenamiento (Site A, no compartidas): %d filas, %d cepas\n",
            length(train_rows), uniqueN(M$strain_id[train_rows])))
print(M[train_rows][, .N, by = mechanism])

fit_and_predict <- function(tr_rows) {
  ytr <- factor(M$mechanism[tr_rows], levels = c("NDM", "KPC"))

  pv <- apply(X[tr_rows, , drop = FALSE], 2,
              function(v) tryCatch(wilcox.test(v ~ ytr)$p.value, error = function(e) 1))
  keep <- which(p.adjust(pv, "fdr") < 0.05)
  if (!length(keep)) keep <- order(pv)[seq_len(min(20, ncol(X)))]
  fit <- randomForest(x = X[tr_rows, keep, drop = FALSE], y = ytr, ntree = 500)

  out <- rbindlist(lapply(pairs_sp, function(s) {
    ra <- which(M$sp == s & M$source_lab == "ALVAREZ")[1]
    rm_ <- which(M$sp == s & M$source_lab == "MALBRAN")[1]
    data.table(
      sp       = s,
      strain_id = M$strain_id[ra],
      truth    = M$mechanism[ra],
      call_A   = as.character(predict(fit, X[ra, keep, drop = FALSE])),
      call_M   = as.character(predict(fit, X[rm_, keep, drop = FALSE]))
    )
  }))
  list(pred = out, n_feat = length(keep))
}

main_fit <- fit_and_predict(train_rows)
pred <- main_fit$pred
cat(sprintf("Features seleccionadas: %d\n\n", main_fit$n_feat))

pred[, concordante := call_A == call_M]
n     <- nrow(pred)
agree <- sum(pred$concordante)
tab   <- table(factor(pred$call_A, c("NDM","KPC")), factor(pred$call_M, c("NDM","KPC")))

cat("=== Tabla pareada (filas = llamada en Site A, columnas = Site M) ===\n")
print(tab)

kappa_of <- function(d) {
  t <- table(factor(d$call_A, c("NDM","KPC")), factor(d$call_M, c("NDM","KPC")))
  nn <- sum(t); po <- sum(diag(t)) / nn
  pe <- sum(rowSums(t) * colSums(t)) / nn^2
  c(agree = po, kappa = (po - pe) / (1 - pe))
}
po <- sum(diag(tab)) / n
pe <- sum(rowSums(tab) * colSums(tab)) / n^2
kappa <- (po - pe) / (1 - pe)

disc_AM <- tab[1, 2]
disc_MA <- tab[2, 1]
mc <- if ((disc_AM + disc_MA) > 0) binom.test(disc_AM, disc_AM + disc_MA, 0.5) else NULL

cat(sprintf("\nConcordancia bruta: %d/%d pares = %.1f%%\n", agree, n, 100 * agree / n))
cat(sprintf("Kappa de Cohen (estimador puntual): %.3f\n", kappa))

pred[, treatment := sub(".*__", "", sp)]
cat("\n=== Por protocolo (n=22 cepas independientes cada uno) ===\n")
for (tt in sort(unique(pred$treatment))) {
  r <- kappa_of(pred[treatment == tt])
  cat(sprintf("  %s: agree=%.3f kappa=%.3f\n", tt, r["agree"], r["kappa"]))
}

byst <- pred[, .(both = all(call_A == call_M)), by = strain_id]
cat(sprintf("\nCepas con ambos protocolos concordantes: %d/%d = %.1f%%\n",
            sum(byst$both), nrow(byst), 100 * mean(byst$both)))

if (!is.null(mc)) {
  cat(sprintf("McNemar exacto (discordantes %d vs %d): p = %.3f\n", disc_AM, disc_MA, mc$p.value))
} else {
  cat("McNemar: sin pares discordantes\n")
}

acc_A <- mean(pred$call_A == pred$truth)
acc_M <- mean(pred$call_M == pred$truth)
cat(sprintf("\nExactitud en Site A: %.3f | en Site M: %.3f (n = %d pares)\n", acc_A, acc_M, n))

tr_strains <- unique(M$strain_id[train_rows])
boot <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  set.seed(1000 + b)
  s  <- sample(tr_strains, length(tr_strains), replace = TRUE)
  rr <- unlist(lapply(s, function(x) train_rows[M$strain_id[train_rows] == x]))
  if (uniqueN(M$mechanism[rr]) < 2) return(NULL)
  p <- tryCatch(fit_and_predict(rr)$pred, error = function(e) NULL)
  if (is.null(p)) return(NULL)
  data.table(iter = b,
             agreement = mean(p$call_A == p$call_M),
             acc_A = mean(p$call_A == p$truth),
             acc_M = mean(p$call_M == p$truth))
}))

set.seed(7)
ev_strains <- unique(pred$strain_id)
cl <- rbindlist(lapply(seq_len(2000), function(i) {
  s <- sample(ev_strains, length(ev_strains), replace = TRUE)
  d <- rbindlist(lapply(s, function(x) pred[strain_id == x]))
  as.data.table(as.list(kappa_of(d)))
}))
ci_agree <- quantile(cl$agree, c(.025, .975))
ci_kappa <- quantile(cl$kappa, c(.025, .975), na.rm = TRUE)
cat("\n=== Cluster bootstrap sobre las cepas evaluadas (2000 iter, IC valido) ===\n")
cat(sprintf("  Concordancia: %.3f (IC95 %.3f a %.3f)\n",
            mean(cl$agree), ci_agree[1], ci_agree[2]))
cat(sprintf("  Kappa:        %.3f (IC95 %.3f a %.3f)\n",
            mean(cl$kappa, na.rm = TRUE), ci_kappa[1], ci_kappa[2]))

cat(sprintf("\n=== Bootstrap del train por cepa (%d iteraciones validas) ===\n", nrow(boot)))
cat(sprintf("Concordancia  media %.3f  IC95 %.3f a %.3f\n",
            mean(boot$agreement), quantile(boot$agreement, .025), quantile(boot$agreement, .975)))
cat(sprintf("Exactitud A   media %.3f  IC95 %.3f a %.3f\n",
            mean(boot$acc_A), quantile(boot$acc_A, .025), quantile(boot$acc_A, .975)))
cat(sprintf("Exactitud M   media %.3f  IC95 %.3f a %.3f\n",
            mean(boot$acc_M), quantile(boot$acc_M, .025), quantile(boot$acc_M, .975)))

save_table(pred, "kpc_ndm_platform_concordance_pairs")
save_table(boot, "kpc_ndm_platform_concordance_bootstrap")
save_table(data.table(
  n_pairs = n, n_strains = uniqueN(pred$strain_id),
  agreement = round(agree / n, 4), kappa = round(kappa, 4),
  agree_lo_cluster = round(ci_agree[1], 4), agree_hi_cluster = round(ci_agree[2], 4),
  kappa_lo_cluster = round(ci_kappa[1], 4), kappa_hi_cluster = round(ci_kappa[2], 4),
  strain_level_both_protocols = round(mean(byst$both), 4),
  mcnemar_p = if (!is.null(mc)) round(mc$p.value, 4) else NA_real_,
  discordant_A_NDM_M_KPC = disc_AM, discordant_A_KPC_M_NDM = disc_MA,
  acc_site_A = round(acc_A, 4), acc_site_M = round(acc_M, 4),
  boot_agree_lo = round(quantile(boot$agreement, .025), 4),
  boot_agree_hi = round(quantile(boot$agreement, .975), 4)
), "kpc_ndm_platform_concordance_summary")

cat("\nGuardado en:", TXT_MODULAR_TAB, "\n")
