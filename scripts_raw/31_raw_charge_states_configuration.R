source("scripts_raw/00_shared_utils_raw.R")

suppressPackageStartupMessages(library(data.table))

OUT <- file.path(dirname(RAW_MODULAR_OUT), "CHARGE_STATES")
dir.create(OUT, showWarnings = FALSE)
PROTON <- 1.007
TOL <- 5e-4
PARES <- list(c(2692.1, 5381.6), c(4154.8, 8309.9), c(5555.2, 11109.3), c(5564.3, 11127.0))
MARCADORES <- c(2416, 2692, 3180, 4154, 3761, 5555, 11109, 11126, 5249, 3348, 5804, 11185)

cliffs <- function(a, b) {
  u <- wilcox.test(a, b, exact = FALSE)$statistic
  as.numeric(2 * u / (length(a) * length(b)) - 1)
}

laser_de <- function(path) {
  h <- readChar(path, 9000, useBytes = TRUE)
  m <- regmatches(h, regexec('value="([0-9.]+)" name="laser energy"', h, useBytes = TRUE))[[1]]
  if (length(m) == 2) as.numeric(m[2]) else NA_real_
}

main <- function() {
  pooled <- readRDS(file.path(RAW_MODULAR_OBJ, "raw_pooled_dataset.rds"))
  meta <- as.data.table(pooled$metadata)
  meta[, laser := vapply(path, laser_de, numeric(1))]
  mz <- as.numeric(colnames(pooled$matrix))
  agg <- aggregate_raw_by_group(pooled$matrix, meta, "strain_treatment")
  m <- as.data.table(agg$metadata); X <- agg$matrix
  nearest <- function(v) { j <- which.min(abs(mz - v)); if (abs(mz[j] - v) <= TOL * v * 2) j else NA_integer_ }
  A <- m$source_lab == "ALVAREZ" & !is.na(m$label_res_sen)
  aa <- which(A); C_A <- suppressWarnings(cor(X[aa, ]))

  filas <- list()
  for (v in MARCADORES) {
    i <- nearest(v); if (is.na(i)) next
    for (z in 2:3) for (dir in c("mas_pesado", "mas_liviano")) {
      target <- if (dir == "mas_pesado") z * mz[i] - (z - 1) * PROTON else (mz[i] + (z - 1) * PROTON) / z
      j <- which(abs(mz - target) <= 2 * TOL * target)
      for (k in j) filas[[length(filas) + 1]] <- data.table(
        marcador = v, mz_bin = round(mz[i], 1), z = z, sentido = dir, par = round(mz[k], 1),
        r = round(C_A[i, k], 3), codeteccion = round(mean(X[aa, i] > 0 & X[aa, k] > 0), 3))
    }
  }
  pares <- rbindlist(filas); setorder(pares, marcador, -r)
  print(pares); fwrite(pares, file.path(OUT, "pares_de_carga_por_marcador.csv"))

  sen <- m[label_res_sen == "SEN"]
  shared <- intersect(sen[source_lab == "ALVAREZ", strain_id], sen[source_lab == "MALBRAN", strain_id])
  cat("Sensibles en ambos sitios:", length(shared), "\n")
  mean_cepa <- function(j, lab) {
    idx <- which(m$label_res_sen == "SEN" & m$source_lab == lab & m$strain_id %in% shared)
    r <- tapply(X[idx, j], m$strain_id[idx], mean); r[shared]
  }
  ef <- list()
  for (p in PARES) {
    i2 <- nearest(p[1]); i1 <- nearest(p[2])
    if (is.na(i2) || is.na(i1)) { cat("sin bin:", p, "\n"); next }
    d <- data.table(z2_new = mean_cepa(i2, "ALVAREZ"), z2_old = mean_cepa(i2, "MALBRAN"),
                    z1_new = mean_cepa(i1, "ALVAREZ"), z1_old = mean_cepa(i1, "MALBRAN"))
    ok <- with(d, (z2_new + z1_new) > 0 & (z2_old + z1_old) > 0)
    frac_new <- with(d[ok], z2_new / (z2_new + z1_new)); frac_old <- with(d[ok], z2_old / (z2_old + z1_old))
    ef[[length(ef) + 1]] <- data.table(
      z2 = p[1], z1 = p[2], cepas = sum(ok),
      z2_old_med = round(median(d$z2_old)), z2_new_med = round(median(d$z2_new)),
      p_z2 = signif(wilcox.test(d$z2_new, d$z2_old, paired = TRUE, exact = FALSE)$p.value, 2),
      z1_old_med = round(median(d$z1_old)), z1_new_med = round(median(d$z1_new)),
      p_z1 = signif(wilcox.test(d$z1_new, d$z1_old, paired = TRUE, exact = FALSE)$p.value, 2),
      frac2_old = round(median(frac_old), 3), frac2_new = round(median(frac_new), 3),
      p_frac = signif(wilcox.test(frac_new, frac_old, paired = TRUE, exact = FALSE)$p.value, 2))
  }
  ef <- rbindlist(ef); print(ef); fwrite(ef, file.path(OUT, "efecto_configuracion_misma_cepa.csv"))

  cl <- list()
  y <- m$label_res_sen[aa]
  for (p in PARES) {
    i2 <- nearest(p[1]); i1 <- nearest(p[2]); if (is.na(i2) || is.na(i1)) next
    fr <- X[aa, i2] / pmax(X[aa, i2] + X[aa, i1], 1)
    cl[[length(cl) + 1]] <- data.table(
      z2 = p[1], z1 = p[2],
      delta_z2 = round(cliffs(X[aa, i2][y == "RES"], X[aa, i2][y == "SEN"]), 3),
      delta_z1 = round(cliffs(X[aa, i1][y == "RES"], X[aa, i1][y == "SEN"]), 3),
      delta_frac2 = round(cliffs(fr[y == "RES"], fr[y == "SEN"]), 3))
  }
  cl <- rbindlist(cl); print(cl); fwrite(cl, file.path(OUT, "clase_sitio_A_z2_z1.csv"))

  las <- list()
  for (p in PARES) {
    i2 <- nearest(p[1]); i1 <- nearest(p[2]); if (is.na(i2) || is.na(i1)) next
    a <- pooled$matrix[, i2]; b <- pooled$matrix[, i1]
    use <- (a + b) > 0 & !is.na(meta$laser)
    fr <- a[use] / (a[use] + b[use]); g <- paste(meta$source_lab[use], meta$label_res_sen[use])
    las[[length(las) + 1]] <- rbindlist(lapply(unique(g), function(k) {
      s <- g == k
      ct <- suppressWarnings(cor.test(fr[s], meta$laser[use][s], method = "spearman", exact = FALSE))
      data.table(z2 = p[1], grupo = k, n = sum(s), rho = round(ct$estimate, 3), p = signif(ct$p.value, 2))
    }))
  }
  las <- rbindlist(las); print(las); fwrite(las, file.path(OUT, "laser_fraccion_2mas.csv"))
  print(meta[, .(n = .N, laser = round(mean(laser, na.rm = TRUE), 3)), by = .(source_lab, label_res_sen)])
  peaklist_fraccion()
}

peaklist_fraccion <- function() {
  e <- new.env(); sys.source("scripts_txt/00_shared_utils.R", envir = e)
  meta <- as.data.table(readRDS(file.path(e$TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
  X <- readRDS(file.path(e$TXT_MODULAR_OBJ, "txt_representations.rds"))$intensity_raw
  ok <- !is.na(meta$label_res_sen); meta <- meta[ok]; X <- X[ok, , drop = FALSE]
  b <- as.numeric(colnames(X)); ventana <- function(c0) rowSums(X[, abs(b - c0) <= 4, drop = FALSE])
  out <- rbindlist(lapply(list(c(2692, 5382), c(4154, 8310), c(5555, 11110)), function(p) {
    z2 <- ventana(p[1]); z1 <- ventana(p[2]); fr <- z2 / pmax(z2 + z1, 1)
    d <- data.table(meta, z2, z1, fr)[(z2 + z1) > 0]
    res <- d[, .(n = .N, z2_med = round(median(z2)), z1_med = round(median(z1)), frac2 = round(median(fr), 3)),
             by = .(source_lab, label_res_sen)][order(source_lab, label_res_sen)]
    m <- d[source_lab == "MALBRAN"]
    sh <- intersect(d[source_lab == "ALVAREZ" & label_res_sen == "RES", strain_id], d[source_lab == "MALBRAN" & label_res_sen == "RES", strain_id])
    a <- d[strain_id %in% sh & source_lab == "ALVAREZ", .(f = mean(fr)), by = strain_id]
    mm <- d[strain_id %in% sh & source_lab == "MALBRAN", .(f = mean(fr)), by = strain_id]
    j <- merge(a, mm, by = "strain_id")
    res[, `:=`(z2 = p[1], z1 = p[2],
               p_clase_sitio_M = signif(wilcox.test(fr ~ label_res_sen, data = m, exact = FALSE)$p.value, 2),
               crkp_compartidas = nrow(j), frac2_crkp_A = round(median(j$f.x), 3), frac2_crkp_M = round(median(j$f.y), 3),
               p_crkp_A_vs_M = signif(wilcox.test(j$f.x, j$f.y, paired = TRUE, exact = FALSE)$p.value, 2))]
  }))
  print(out); fwrite(out, file.path(OUT, "peaklist_fraccion_2mas_por_sitio.csv"))
}

main()
