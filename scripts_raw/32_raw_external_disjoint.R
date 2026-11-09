source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
})

N_ITER <- 100
OUT <- file.path(dirname(RAW_MODULAR_OUT), "SITE_M")
dir.create(OUT, showWarnings = FALSE)

fit_rf <- function(X, y, balanced) {
  ytr <- factor(y)
  if (!balanced) return(randomForest(x = X, y = ytr, ntree = 500))
  n <- min(table(ytr))
  randomForest(x = X, y = ytr, ntree = 500, strata = ytr, sampsize = c(n, n))
}

main <- function() {
  set.seed(42)
  pooled <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_pooled_dataset.rds"))
  agg <- aggregate_raw_by_group(pooled$matrix, as.data.table(pooled$metadata), "strain_treatment")
  m <- as.data.table(agg$metadata); X <- agg$matrix
  tr <- which(m$source_lab == "ALVAREZ" & !is.na(m$label_res_sen))
  te <- which(m$source_lab == "MALBRAN" & !is.na(m$label_res_sen))
  sen_a <- unique(m[tr][label_res_sen == "SEN", strain_id])
  sen_m <- unique(m[te][label_res_sen == "SEN", strain_id])
  shared <- intersect(sen_a, sen_m); solo_m <- setdiff(sen_m, sen_a)
  cat("SEN A:", length(sen_a), "M:", length(sen_m), "compartidas:", length(shared), "solo M:", length(solo_m), "\n")

  filas <- list()
  for (bal in c(FALSE, TRUE)) {
    fit <- fit_rf(X[tr, ], m$label_res_sen[tr], bal)
    pred <- as.character(predict(fit, X[te, ]))
    for (g in c("compartidas", "solo_M", "todas")) {
      idx <- switch(g, compartidas = m$strain_id[te] %in% shared, solo_M = m$strain_id[te] %in% solo_m, todas = rep(TRUE, length(te)))
      k <- sum(pred[idx] == "SEN"); n <- sum(idx); ic <- binom.test(k, n)$conf.int
      filas[[length(filas) + 1]] <- data.table(balanceado = bal, conjunto = g, cepas = uniqueN(m$strain_id[te][idx]),
        obs = n, llamadas_sen = k, especificidad = round(k / n, 3), lo = round(ic[1], 3), hi = round(ic[2], 3))
    }
  }
  r1 <- rbindlist(filas); print(r1); fwrite(r1, file.path(OUT, "raw_especificidad_externa_por_conjunto.csv"))

  it <- rbindlist(lapply(seq_len(N_ITER), function(i) rbindlist(lapply(c(TRUE, FALSE), function(bal) {
    half <- sample(shared, floor(length(shared) / 2))
    other <- setdiff(sen_m, half)
    trn <- tr[m$label_res_sen[tr] == "RES" | m$strain_id[tr] %in% half]
    tst <- te[m$strain_id[te] %in% other]
    fit <- fit_rf(X[trn, ], m$label_res_sen[trn], bal)
    data.table(iter = i, balanceado = bal, cepas_sen_train = length(half), cepas_test = length(other),
               spec = mean(as.character(predict(fit, X[tst, ])) == "SEN"))
  }))))
  r2 <- it[, .(iter = .N, cepas_sen_train = first(cepas_sen_train), cepas_test = first(cepas_test),
               spec = round(mean(spec), 3), sd = round(sd(spec), 3),
               lo = round(quantile(spec, .025), 3), hi = round(quantile(spec, .975), 3)), by = balanceado]
  print(r2); fwrite(r2, file.path(OUT, "raw_especificidad_disjunta.csv"))
}

main()
