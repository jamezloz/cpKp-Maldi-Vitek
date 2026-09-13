source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({ library(data.table) })

RAW_TAB <- "resultados/RAW_modular/tables"
RAW_OBJ <- "resultados/RAW_modular/objects"

STAB_THRESH  <- 0.8
LARGE_DELTA  <- 0.474

MATCH_DA     <- 10
EXACT_DA     <- 0.5
N_BOOT_CI    <- 1000

PROPUESTOS_P1 <- list(res_sen = c(2416, 2692, 3180, 4154),
                   kpc_ndm = c(3761, 5555, 11110, 11126))

read_if <- function(p) if (file.exists(p)) fread(p) else NULL

cliff_fast <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  n1 <- length(a); n2 <- length(b)
  if (!n1 || !n2) return(NA_real_)
  r <- rank(c(a, b))
  auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
  2 * auc - 1
}

cliff_ci <- function(va, vb, sa, sb, n_boot = N_BOOT_CI) {
  point <- cliff_fast(va, vb)
  if (is.na(point)) return(list(delta = NA_real_, lo = NA_real_, hi = NA_real_))
  ua <- unique(sa); ub <- unique(sb)
  ia <- split(seq_along(va), sa); ib <- split(seq_along(vb), sb)
  boots <- numeric(n_boot)
  for (k in seq_len(n_boot)) {
    ra <- unlist(ia[sample(ua, replace = TRUE)], use.names = FALSE)
    rb <- unlist(ib[sample(ub, replace = TRUE)], use.names = FALSE)
    boots[k] <- cliff_fast(va[ra], vb[rb])
  }
  q <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  list(delta = point, lo = unname(q[1]), hi = unname(q[2]))
}

nearest_idx <- function(mz_target, mz_list, tol_da = MATCH_DA) {
  d <- abs(mz_list - mz_target)
  i <- which.min(d)
  if (length(i) && d[i] <= tol_da) i else NA_integer_
}

nearest_val <- function(mz_target, mz_list, values, tol_da = MATCH_DA) {
  i <- nearest_idx(mz_target, mz_list, tol_da)
  if (is.na(i)) NA_real_ else values[i]
}

txt_mat_raw  <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_matrix_aggregated.rds"))
txt_meta_raw <- as.data.table(readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds")))
{
  key <- paste(txt_meta_raw$strain_id, txt_meta_raw$treatment, sep = "__")
  idx <- split(seq_len(nrow(txt_meta_raw)), key)
  txt_mat <- t(vapply(idx, function(i) colMeans(txt_mat_raw[i, , drop = FALSE]),
                      numeric(ncol(txt_mat_raw))))
  txt_meta <- rbindlist(lapply(idx, function(i) txt_meta_raw[i][1,
                          .(strain_id, treatment, mechanism, mechanism_raw, label_res_sen)]))
}
txt_col_mz <- as.numeric(colnames(txt_mat))

raw_ds   <- readRDS(file.path(RAW_OBJ, "raw_peak_strain_treatment_dataset.rds"))
raw_mat  <- raw_ds$matrix
raw_meta <- as.data.table(raw_ds$metadata)
raw_col_mz <- as.numeric(colnames(raw_mat))

task_groups <- list(
  res_sen = list(txt = function(m) ifelse(m$label_res_sen == "RES", "RES", "SEN"),
                 raw = function(m) ifelse(m$label_res_sen == "RES", "RES", "SEN"),
                 a = "RES", b = "SEN"),
  kpc_ndm = list(txt = function(m) m$mechanism,
                 raw = function(m) m$mechanism,
                 a = "KPC", b = "NDM")
)

build_task <- function(task) {
  tg <- task_groups[[task]]

  stab_txt <- read_if(file.path(TXT_MODULAR_TAB, sprintf("stability_selection_%s.csv", task)))
  stab_raw <- read_if(file.path(RAW_TAB, "improved_peaklist",
                                sprintf("bootstrap_stability_%s.csv", task)))
  raw_uni  <- read_if(file.path(RAW_TAB, sprintf("raw_classic_final_biomarkers_%s.csv", task)))
  shap_txt <- read_if(file.path(TXT_MODULAR_TAB, sprintf("shap_importance_%s.csv", task)))
  persite  <- read_if(file.path(TXT_MODULAR_TAB, "per_hospital_stability.csv"))
  if (is.null(stab_txt) || is.null(stab_raw)) {
    warning("Faltan estabilidades para ", task, call. = FALSE, immediate. = TRUE); return(NULL)
  }
  stab_txt[, feature := as.numeric(feature)]
  stab_raw[, mz := as.numeric(mz)]

  shap_mz <- if (!is.null(shap_txt)) as.numeric(shap_txt$feature) else numeric(0)
  shap_rk <- if (!is.null(shap_txt)) as.numeric(seq_len(nrow(shap_txt))) else numeric(0)

  shap_raw <- read_if(file.path(RAW_TAB, sprintf("shap_importance_%s_raw.csv", task)))
  shap_raw_mz <- if (!is.null(shap_raw)) as.numeric(shap_raw$feature) else numeric(0)
  shap_raw_rk <- if (!is.null(shap_raw)) as.numeric(seq_len(nrow(shap_raw))) else numeric(0)

  cand_src <- rbindlist(list(
    stab_txt[selection_freq >= STAB_THRESH, .(mz = feature, w = selection_freq)],
    stab_raw[stability     >= STAB_THRESH, .(mz = mz,      w = stability)]
  ))
  setorder(cand_src, mz)

  pubs <- sort(PROPUESTOS_P1[[task]])
  far  <- cand_src[vapply(mz, function(m) all(abs(pubs - m) > MATCH_DA),
                          logical(1))]
  reps_far <- numeric(0)
  if (nrow(far)) {
    gf <- cumsum(c(1L, as.integer(diff(far$mz) > MATCH_DA)))
    reps_far <- vapply(split(seq_len(nrow(far)), gf),
                       function(ii) far$mz[ii][which.max(far$w[ii])], numeric(1))
  }
  cand <- sort(unique(c(pubs, reps_far)))

  rows <- lapply(cand, function(mz) {
    jt <- nearest_idx(mz, txt_col_mz)
    jr <- nearest_idx(mz, raw_col_mz)

    st <- if (!is.na(jt)) nearest_val(txt_col_mz[jt], stab_txt$feature, stab_txt$selection_freq, EXACT_DA) else NA_real_
    sr <- if (!is.na(jr)) nearest_val(raw_col_mz[jr], stab_raw$mz,      stab_raw$stability,      MATCH_DA) else NA_real_

    ct <- list(delta = NA_real_, lo = NA_real_, hi = NA_real_)
    if (!is.na(jt)) {
      gt <- tg$txt(txt_meta); kt <- gt %in% c(tg$a, tg$b)
      ct <- cliff_ci(txt_mat[kt & gt == tg$a, jt], txt_mat[kt & gt == tg$b, jt],
                     txt_meta$strain_id[kt & gt == tg$a], txt_meta$strain_id[kt & gt == tg$b])
    }
    cr <- list(delta = NA_real_, lo = NA_real_, hi = NA_real_)
    if (!is.na(jr)) {
      gr <- tg$raw(raw_meta); kr <- gr %in% c(tg$a, tg$b)
      cr <- cliff_ci(raw_mat[kr & gr == tg$a, jr], raw_mat[kr & gr == tg$b, jr],
                     raw_meta$strain_id[kr & gr == tg$a], raw_meta$strain_id[kr & gr == tg$b])
    }

    data.table(
      tarea = task, mz = round(mz, 1),
      estab_TXT = st, estab_RAW = sr,
      cliff_TXT = ct$delta, cliff_TXT_lo = ct$lo, cliff_TXT_hi = ct$hi,
      cliff_RAW = cr$delta, cliff_RAW_lo = cr$lo, cliff_RAW_hi = cr$hi,
      shap_rank_TXT = nearest_val(mz, shap_mz, shap_rk),
      shap_rank_RAW = nearest_val(mz, shap_raw_mz, shap_raw_rk),
      propuesto_p1 = any(abs(PROPUESTOS_P1[[task]] - mz) <= MATCH_DA)
    )
  })
  dt <- rbindlist(rows)

  dt[, estable_TXT := !is.na(estab_TXT) & estab_TXT >= STAB_THRESH]
  dt[, estable_RAW := !is.na(estab_RAW) & estab_RAW >= STAB_THRESH]
  dt[, dir_consistente := !is.na(cliff_TXT) & !is.na(cliff_RAW) &
                          sign(cliff_TXT) == sign(cliff_RAW)]

  if (!is.null(persite)) {
    a <- persite[sitio == "ALVAREZ" & tarea == task]
    m <- persite[sitio == "MALBRAN" & tarea == task]
    match_site <- function(d) vapply(dt$mz, function(mz) {
      if (!nrow(d)) return(NA_real_)
      i <- which.min(abs(as.numeric(d$feature) - mz))
      if (length(i) && abs(as.numeric(d$feature)[i] - mz) <= MATCH_DA) d$selection_freq[i] else NA_real_
    }, numeric(1))
    dt[, estab_Alvarez := match_site(a)]
    dt[, estab_Malbran := match_site(m)]
  } else {
    warning("per_hospital_stability ausente: sin columnas de sitio", call. = FALSE, immediate. = TRUE)
  }

  dt[, efecto_grande_RAW := !is.na(cliff_RAW) & abs(cliff_RAW) >= LARGE_DELTA]

  dt[, ic_TXT_solido := !is.na(cliff_TXT_lo) & !is.na(cliff_TXT_hi) & cliff_TXT_lo * cliff_TXT_hi > 0]
  dt[, ic_RAW_solido := !is.na(cliff_RAW_lo) & !is.na(cliff_RAW_hi) & cliff_RAW_lo * cliff_RAW_hi > 0]

  dt[, nivel := fifelse(estable_TXT & estable_RAW & dir_consistente, "1_ambas_plataformas",
                fifelse(estable_RAW, "2_solo_RAW",
                fifelse(estable_TXT, "3_solo_TXT", "4_no_estable")))]

  dt[, panel_principal := nivel == "1_ambas_plataformas" & efecto_grande_RAW &
                          ic_TXT_solido & ic_RAW_solido]

  dt[, efecto_grande_TXT := !is.na(cliff_TXT) & abs(cliff_TXT) >= LARGE_DELTA]
  dt[, fuerte_peaklist := nivel == "1_ambas_plataformas" & efecto_grande_TXT &
                          ic_TXT_solido & !panel_principal]

  dt[, abs_cliff := pmax(abs(cliff_TXT), abs(cliff_RAW), na.rm = TRUE)]
  setorder(dt, nivel, -abs_cliff, na.last = TRUE)
  dt[]
}

main <- function() {
  out <- rbindlist(lapply(c("res_sen", "kpc_ndm"), build_task), fill = TRUE)
  save_table(out, "biomarker_panel_consolidado")

  for (tk in unique(out$tarea)) {
    d <- out[tarea == tk]
    cat("\n==============", tk, "==============\n")
    cat("reparto por nivel:\n"); print(d[, .N, by = nivel][order(nivel)])
    cat("panel principal:", d[panel_principal == TRUE, .N], "masas\n")
    cc <- c("mz", "estab_TXT", "estab_RAW", "cliff_TXT", "cliff_RAW",
            "dir_consistente", "propuesto_p1", "nivel")
    cat("\npropuestas en paper 1:\n"); print(d[propuesto_p1 == TRUE, ..cc])
  }
  cat("\nTabla en:", file.path(TXT_MODULAR_TAB, "biomarker_panel_consolidado.csv"), "\n")
}

set.seed(42)
main()
