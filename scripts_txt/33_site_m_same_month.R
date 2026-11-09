source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
})

N_REP <- 200
OUT <- file.path(dirname(TXT_MODULAR_OUT), "SITE_M")
dir.create(OUT, showWarnings = FALSE)
mes_de_cabecera <- function(path) {
  h <- readLines(path, n = 5, warn = FALSE)
  m <- regmatches(h, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", h))
  if (length(m)) substr(m[1], 1, 7) else NA_character_
}

predecir <- function(X, y, tr, te) {
  ytr <- factor(y[tr]); n <- min(table(ytr))
  fit <- randomForest(x = X[tr, , drop = FALSE], y = ytr, ntree = 500, strata = ytr, sampsize = c(n, n))
  predict(fit, X[te, , drop = FALSE], type = "prob")[, "RES"]
}

resumir <- function(p, y) {
  sens <- mean(p[y == "RES"] > 0.5); spec <- mean(p[y == "SEN"] <= 0.5)
  c(sens = sens, spec = spec, BA = (sens + spec) / 2)
}

main <- function() {
  set.seed(42)
  meta <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_raw.rds")))
  X <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_raw.rds"))
  stopifnot(nrow(meta) == nrow(X))
  meta[, fila := .I]
  meta <- meta[source_lab == "MALBRAN" & !is.na(label_res_sen)]
  meta[, mes := vapply(path, mes_de_cabecera, character(1))]
  cat("Espectros del sitio M con mes:", nrow(meta), "\n")

  g <- meta[, .(filas = list(fila)), by = .(strain_id, treatment, mes, label_res_sen)]
  A <- t(vapply(g$filas, function(i) colMeans(X[i, , drop = FALSE]), numeric(ncol(X))))
  y <- g$label_res_sen
  print(g[, .(obs = .N, cepas = uniqueN(strain_id)), by = .(mes, label_res_sen)][order(mes)])

  nov <- unique(g[mes == "2024-11", strain_id])
  te <- which(g$mes == "2024-11"); tr <- which(!g$strain_id %in% nov)
  k_res <- uniqueN(g$strain_id[te][y[te] == "RES"]); k_sen <- uniqueN(g$strain_id[te][y[te] == "SEN"])
  obs <- rowMeans(vapply(seq_len(N_REP), function(i) resumir(predecir(A, y, tr, te), y[te]), numeric(3)))

  cepas <- unique(g[, .(strain_id, label_res_sen)])
  nulo <- t(vapply(seq_len(N_REP), function(i) {
    s <- c(sample(cepas[label_res_sen == "RES", strain_id], k_res), sample(cepas[label_res_sen == "SEN", strain_id], k_sen))
    te0 <- which(g$strain_id %in% s); tr0 <- which(!g$strain_id %in% s)
    resumir(predecir(A, y, tr0, te0), y[te0])
  }, numeric(3)))

  res <- data.table(diseno = c("noviembre_2024", "azar_mismo_tamano"),
                    cepas_test_RES = k_res, cepas_test_SEN = k_sen,
                    cepas_train_RES = uniqueN(g$strain_id[tr][y[tr] == "RES"]),
                    sens = round(c(obs["sens"], mean(nulo[, "sens"])), 3),
                    spec = round(c(obs["spec"], mean(nulo[, "spec"])), 3),
                    BA = round(c(obs["BA"], mean(nulo[, "BA"])), 3),
                    sens_p025 = c(NA, round(quantile(nulo[, "sens"], .025), 3)),
                    sens_p975 = c(NA, round(quantile(nulo[, "sens"], .975), 3)),
                    frac_azar_peor_sens = c(NA, round(mean(nulo[, "sens"] <= obs["sens"]), 3)))
  print(res); fwrite(res, file.path(OUT, "control_mismo_mes.csv"))
}

main()
