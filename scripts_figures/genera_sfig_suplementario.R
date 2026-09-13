source("scripts_figures/00_fig_utils.R")

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
WHICH  <- if (length(args) >= 1) as.integer(args[1]) else 0L
do_fig <- function(n) WHICH == 0L || WHICH == n

COL_PHENO <- c(CRKp = "#C0392B", CSKp = "#2980B9")
TASK_LAB  <- c(res_sen = "CRKp/CSKp", kpc_ndm = "KPC/NDM")

fig_s4 <- function() {
  cat("== Figure S4: bootstrap stability curves ==\n")
  load_txt <- function(f) {
    d <- fread(f); setnames(d, c("feature", "selection_freq"), c("mz", "stability"))
    d[, mz := as.numeric(mz)][]
  }
  load_raw <- function(f) fread(f)[, .(mz = as.numeric(mz), stability)]

  dt <- rbindlist(list(
    load_txt(file.path(TXT_TAB, "stability_selection_res_sen.csv"))[, `:=`(panel = "A  peak-list — CRKp/CSKp")],
    load_txt(file.path(TXT_TAB, "stability_selection_kpc_ndm.csv"))[, `:=`(panel = "B  peak-list — KPC/NDM")],
    load_raw(file.path(RAW_TAB, "improved_peaklist/bootstrap_stability_res_sen.csv"))[, `:=`(panel = "C  RAW dataset — CRKp/CSKp")],
    load_raw(file.path(RAW_TAB, "improved_peaklist/bootstrap_stability_kpc_ndm.csv"))[, `:=`(panel = "D  RAW dataset — KPC/NDM")]
  ), use.names = TRUE)
  dt[, stable := stability >= 0.9999]

  p <- ggplot(dt, aes(mz, stability)) +
    geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    geom_point(aes(colour = stable), size = 1.1, alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "Stability = 1.00", `FALSE` = "< 1.00"),
                        name = NULL) +
    facet_wrap(~ panel, scales = "free_x", ncol = 2) +
    labs(x = expression(italic("m/z")~"(Da)"), y = "Bootstrap selection frequency") +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    theme_paper(12) + theme(legend.position = "top")

  save_fig(p, "sfig4_stability_curves_suplementario", width = 26, height = 20)
}

fig_s6 <- function() {
  cat("== Figure S6: BA/AUC distributions ==\n")
  raw <- fread(file.path(RAW_TAB, "raw_extended_classifiers_strain_treatment.csv"))
  txt <- fread(file.path(TXT_TAB, "ml_extended_classifiers_iterations.csv"))
  keep <- c("model", "task", "balanced_accuracy", "auc")
  raw <- raw[, ..keep][, dataset := "RAW dataset"]
  txt <- txt[, ..keep][, dataset := "peak-list"]
  d <- rbindlist(list(txt, raw), use.names = TRUE)

  MODEL_LAB <- c(rf = "RF", xgb = "XGB", svm = "SVM", knn = "KNN", dt = "DT",
                 extra_trees = "Extra Trees", adaboost = "AdaBoost",
                 logistic = "Logistic", naive_bayes = "Naïve Bayes")
  d[, model := MODEL_LAB[model]]
  d[, row := factor(paste0(dataset, " — ", TASK_LAB[task]),
                    levels = c("peak-list — CRKp/CSKp", "peak-list — KPC/NDM",
                               "RAW dataset — CRKp/CSKp", "RAW dataset — KPC/NDM"))]
  long <- melt(d, id.vars = c("model", "row"),
               measure.vars = c("balanced_accuracy", "auc"),
               variable.name = "metric", value.name = "value")
  long[, metric := factor(metric, levels = c("balanced_accuracy", "auc"),
                          labels = c("Balanced accuracy", "AUC"))]

  ord <- d[, .(m = median(balanced_accuracy, na.rm = TRUE)), by = .(row, model)]

  mk <- function(rw) {
    dd <- long[row == rw]
    om <- ord[row == rw][order(m), model]
    dd[, model := factor(model, levels = om)]
    ggplot(dd, aes(model, value)) +
      geom_violin(aes(fill = metric), colour = NA, alpha = 0.45, scale = "width") +
      geom_boxplot(width = 0.18, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
      facet_wrap(~ metric, nrow = 1) +
      coord_flip() +
      scale_fill_manual(values = c("Balanced accuracy" = "#2980B9", "AUC" = "#27AE60"), guide = "none") +
      labs(title = rw, x = NULL, y = NULL) +
      theme_paper(10) + theme(plot.title = element_text(size = 11, face = "bold"))
  }
  rows <- levels(long$row)
  p <- (mk(rows[1]) / mk(rows[2]) / mk(rows[3]) / mk(rows[4]))
  save_fig(p, "sfig6_ba_auc_distributions_suplementario", width = 24, height = 32)
}

fig_s8 <- function() {
  cat("== Figure S8: decision curve analysis ==\n")
  load_dca <- function(f, panel) {
    d <- fread(f)
    melt(d, id.vars = "threshold",
         measure.vars = c("net_benefit_rf", "net_benefit_treat_all", "net_benefit_treat_none"),
         variable.name = "strategy", value.name = "nb")[, panel := panel][]
  }
  dt <- rbindlist(list(
    load_dca(file.path(TXT_TAB, "dca_res_sen_txt.csv"), "A  peak-list — CRKp/CSKp"),
    load_dca(file.path(TXT_TAB, "dca_kpc_ndm_txt.csv"), "B  peak-list — KPC/NDM"),
    load_dca(file.path(RAW_TAB, "dca_res_sen_raw.csv"), "C  RAW dataset — CRKp/CSKp"),
    load_dca(file.path(RAW_TAB, "dca_kpc_ndm_raw.csv"), "D  RAW dataset — KPC/NDM")
  ))
  dt[, strategy := factor(strategy,
        levels = c("net_benefit_rf", "net_benefit_treat_all", "net_benefit_treat_none"),
        labels = c("RF classifier", "Treat all", "Treat none"))]

  p <- ggplot(dt, aes(threshold, nb, colour = strategy, linetype = strategy)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = c("RF classifier" = "#2980B9",
                                   "Treat all" = "grey45", "Treat none" = "grey10"), name = NULL) +
    scale_linetype_manual(values = c("RF classifier" = "solid",
                                     "Treat all" = "dashed", "Treat none" = "dotted"), name = NULL) +
    coord_cartesian(ylim = c(-0.05, NA)) +
    labs(x = "Threshold probability", y = "Net benefit") +
    theme_paper(12) + theme(legend.position = "top")

  save_fig(p, "sfig8_dca_suplementario", width = 26, height = 20)
}

get_dataset <- function(which) {
  if (which == "peak-list") {
    meta <- as.data.table(readRDS(file.path(TXT_OBJ, "txt_metadata_aggregated.rds")))
    d_cos <- readRDS(file.path(TXT_OBJ, "distance_intensity_raw_cosine.rds"))
    d_bray <- readRDS(file.path(TXT_OBJ, "distance_intensity_tic_bray.rds"))
  } else {
    ds <- readRDS(file.path(RAW_OBJ, "raw_peak_strain_treatment_dataset.rds"))
    meta <- as.data.table(ds$metadata)
    dl <- readRDS(file.path(RAW_OBJ, "raw_classic_distances_strain_treatment.rds"))
    d_cos <- dl$peak_intensity_raw_cosine
    d_bray <- dl$peak_intensity_tic_bray
  }
  meta[, phenotype := ifelse(label_res_sen == "RES", "CRKp", "CSKp")]
  meta[, mech := mechanism]
  list(meta = meta, cos = d_cos, bray = d_bray)
}

SHAPE_MECH <- c(KPC = 16, NDM = 17, "KPC-NDM" = 18, SEN = 15)

fig_s1 <- function() {
  cat("== Figure S1: ordination (NMDS/CAP/dbRDA) ==\n")
  suppressPackageStartupMessages(library(vegan))
  panel <- function(scores_dt, meta, xlab, ylab, tag, subtitle) {
    dd <- cbind(scores_dt, meta[, .(phenotype, mech)])
    ggplot(dd, aes(Dim1, Dim2, colour = phenotype, shape = mech)) +
      stat_ellipse(aes(group = phenotype), type = "norm", level = 0.95,
                   linewidth = 0.4, show.legend = FALSE) +
      geom_point(size = 1.5, alpha = 0.8) +
      scale_colour_manual(values = COL_PHENO, name = "Phenotype") +
      scale_shape_manual(values = SHAPE_MECH, name = "Carbapenemase") +
      labs(x = xlab, y = ylab, tag = tag, subtitle = subtitle) +
      theme_paper(11) + theme(plot.subtitle = element_text(size = 10, face = "bold"))
  }
  build <- function(which, tags) {
    ds <- get_dataset(which); meta <- ds$meta
    set.seed(42)
    nm <- metaMDS(ds$bray, k = 2, trace = 0, try = 20, trymax = 50)
    p_nmds <- panel(data.table(Dim1 = nm$points[,1], Dim2 = nm$points[,2]), meta,
                    "NMDS1", "NMDS2", tags[1],
                    sprintf("%s — NMDS (stress = %.3f)", which, nm$stress))
    cap <- capscale(ds$cos ~ phenotype, data = meta)
    cs  <- scores(cap, display = "sites", choices = 1:2)
    p_cap <- panel(data.table(Dim1 = cs[,1], Dim2 = cs[,2]), meta,
                   "CAP1", "MDS1", tags[2], sprintf("%s — CAP (~ phenotype)", which))
    db <- dbrda(ds$cos ~ phenotype + mech, data = meta)
    dbs <- scores(db, display = "sites", choices = 1:2)
    p_db <- panel(data.table(Dim1 = dbs[,1], Dim2 = dbs[,2]), meta,
                  "dbRDA1", "dbRDA2", tags[3], sprintf("%s — dbRDA (~ phenotype + type)", which))
    list(p_nmds, p_cap, p_db)
  }
  pl <- build("peak-list", c("A", "C", "E"))
  pr <- build("RAW dataset", c("B", "D", "F"))
  p <- (pl[[1]] | pr[[1]]) / (pl[[2]] | pr[[2]]) / (pl[[3]] | pr[[3]]) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_fig(p, "sfig1_ordination_suplementario", width = 26, height = 32)
}

fig_s3 <- function() {
  cat("== Figure S3: factorial PERMANOVA ==\n")
  suppressPackageStartupMessages(library(vegan))
  meta <- as.data.table(readRDS(file.path(TXT_OBJ, "txt_metadata_aggregated.rds")))
  reps <- c(raw_cosine = "distance_intensity_raw_cosine.rds",
            tic_bray   = "distance_intensity_tic_bray.rds",
            max_euclidean = "distance_intensity_max_euclidean.rds",
            clr_euclidean = "distance_intensity_clr_euclidean.rds",
            binary_jaccard = "distance_binary_jaccard.rds")
  rows <- list()
  for (nm in names(reps)) {
    d <- readRDS(file.path(TXT_OBJ, reps[nm]))
    set.seed(42)
    a <- adonis2(d ~ treatment + label_res_sen + source_lab, data = meta,
                 by = "margin", permutations = 999)
    a <- as.data.table(a, keep.rownames = "term")
    a <- a[term %in% c("treatment", "label_res_sen", "source_lab")]
    rows[[nm]] <- a[, .(representation = nm, term, R2, p = `Pr(>F)`)]
  }
  dt <- rbindlist(rows)
  dt[, term := factor(term, levels = c("label_res_sen", "treatment", "source_lab"),
                      labels = c("Resistance (CRKp/CSKp)", "Preparation (COL/EXT)", "Acquisition site"))]
  dt[, representation := factor(representation, levels = names(reps))]
  dt[, sig := ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))]

  p <- ggplot(dt, aes(representation, R2, fill = term)) +
    geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
    geom_text(aes(label = sig), position = position_dodge(0.8), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Resistance (CRKp/CSKp)" = "#C0392B",
                                 "Preparation (COL/EXT)" = "#E67E22",
                                 "Acquisition site" = "#2980B9"), name = "Factor") +
    labs(x = "Spectral representation (peak-list)", y = expression(PERMANOVA~R^2~"(marginal)")) +
    theme_paper(12) + theme(legend.position = "top",
                            axis.text.x = element_text(angle = 20, hjust = 1))
  save_fig(p, "sfig3_factorial_permanova_suplementario", width = 24, height = 16)
}

fig_s7 <- function() {
  cat("== Figure S7: SHAP dependence plots for stable biomarkers (peak-list) ==\n")

  bm <- list(
    res_sen = data.table(mz = c(2416, 2692, 3180, 4154),
                         gene = c("rplJ (uL10)", "Ead/Ea22-like", "rpmJ/bL36", "Type IV pilus")),
    kpc_ndm = data.table(mz = c(3762, 5555, 11110, 11126),
                         gene = c("uncharacterized", "uncharacterized", "ureA / p019", "YsaB lipoprotein")))
  pos_lab <- c(res_sen = "CRKp", kpc_ndm = "KPC")
  build <- function(task, tasklab) {
    sh <- fread(file.path(TXT_TAB, paste0("shap_matrix_", task, ".csv")), header = TRUE)
    fv <- fread(file.path(TXT_TAB, paste0("shap_feature_values_", task, ".csv")), header = TRUE)
    cols <- suppressWarnings(as.numeric(names(sh)))
    rows <- list()
    for (i in seq_len(nrow(bm[[task]]))) {
      mz <- bm[[task]]$mz[i]; j <- which.min(abs(cols - mz)); cn <- names(sh)[j]
      rows[[i]] <- data.table(
        intensity = as.numeric(fv[[cn]]), shap = as.numeric(sh[[cn]]),
        panel = sprintf("%s  ~%s Da\n%s", tasklab, cn, bm[[task]]$gene[i]))
    }
    rbindlist(rows)
  }
  dd <- rbind(build("res_sen", "CRKp/CSKp"), build("kpc_ndm", "KPC/NDM"))
  dd[, panel := factor(panel, levels = unique(panel))]
  p <- ggplot(dd, aes(intensity, shap)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(aes(colour = shap > 0), size = 1, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.6, span = 1) +
    scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2980B9"),
                        labels = c(`TRUE` = "pushes toward positive class", `FALSE` = "toward negative"),
                        name = NULL) +
    facet_wrap(~ panel, scales = "free", ncol = 4) +
    scale_x_continuous(n.breaks = 3, labels = function(x) format(x, scientific = TRUE, digits = 1)) +
    labs(x = "Peak intensity (a.u.)", y = "SHAP value") +
    theme_paper(10) + theme(legend.position = "top",
                            strip.text = element_text(size = 8.5, face = "bold"),
                            axis.text.x = element_text(size = 8))
  save_fig(p, "sfig7_shap_dependence_suplementario", width = 30, height = 16)
}

fig_s2 <- function() {
  cat("== Figure S2: dendrograms (NJ / UPGMA) ==\n")
  suppressPackageStartupMessages(library(ape))
  grp_col <- function(meta) {
    g <- ifelse(meta$mechanism == "SEN", "CSKp",
         ifelse(meta$mechanism == "KPC", "CRKp-KPC",
         ifelse(meta$mechanism == "NDM", "CRKp-NDM", "CRKp-KPC+NDM")))
    pal <- c("CRKp-KPC" = "#C0392B", "CRKp-NDM" = "#2980B9",
             "CSKp" = "#7F8C8D", "CRKp-KPC+NDM" = "#8E44AD")
    list(g = g, col = pal[g], pal = pal)
  }
  path <- file.path(FIG_OUT, "sfig2_dendrograms_suplementario.png")
  png(path, width = 26, height = 26, units = "cm", res = 300, bg = "white")
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1), xpd = NA)
  for (which in c("peak-list", "RAW dataset")) {
    ds <- get_dataset(which); meta <- ds$meta; d <- ds$cos
    gc <- grp_col(meta)

    tr <- ape::nj(d); tr$tip.label <- rep("", length(tr$tip.label))
    ord <- match(attr(d, "Labels"), meta$sample)
    plot(tr, type = "fan", show.tip.label = FALSE, edge.width = 0.4,
         main = sprintf("%s — Neighbour-joining", which))
    tiplabels(pch = 16, col = gc$col[ord], cex = 0.5)

    hc <- hclust(d, method = "average")
    coph <- cor(cophenetic(hc), d)
    phy <- as.phylo(hc); phy$tip.label <- rep("", length(phy$tip.label))
    plot(phy, type = "fan", show.tip.label = FALSE, edge.width = 0.4,
         main = sprintf("%s — UPGMA (cophenetic r = %.2f)", which, coph))
    tiplabels(pch = 16, col = gc$col[match(hc$labels, meta$sample)], cex = 0.5)
  }

  pal <- grp_col(get_dataset("peak-list")$meta)$pal
  legend("bottomright", inset = c(-0.05, -0.08), legend = names(pal), col = pal,
         pch = 16, bty = "n", cex = 0.9, ncol = 2)
  dev.off()
  cat("Saved:", path, "\n")
}

fig_s5 <- function(n_iter = 25) {
  cat("== Figure S5: ROC curves (recompute, grouped CV) ==\n")
  suppressPackageStartupMessages({
    library(e1071); library(randomForest); library(ada); library(xgboost); library(pROC); library(glmnet)
  })
  fit_prob <- function(m, Xtr, ytr, Xte, pos) {
    sc <- scale(Xtr); ctr <- attr(sc,"scaled:center"); scl <- attr(sc,"scaled:scale"); scl[scl==0] <- 1
    Xtr_s <- scale(Xtr, ctr, scl); Xte_s <- scale(Xte, ctr, scl)
    if (m == "svm") {
      f <- svm(Xtr_s, ytr, kernel="radial", probability=TRUE)
      attr(predict(f, Xte_s, probability=TRUE), "probabilities")[, pos]
    } else if (m == "rf") {
      f <- randomForest(Xtr, ytr, ntree=500); predict(f, Xte, type="prob")[, pos]
    } else if (m == "extra_trees") {
      f <- randomForest(Xtr, ytr, ntree=500, mtry=ncol(Xtr), replace=FALSE); predict(f, Xte, type="prob")[, pos]
    } else if (m == "adaboost") {
      Xa <- as.data.frame(Xtr); colnames(Xa) <- paste0("f", seq_len(ncol(Xa)))
      Xb <- as.data.frame(Xte); colnames(Xb) <- paste0("f", seq_len(ncol(Xb)))
      f <- ada(Xa, ytr, iter=60, loss="exponential")
      predict(f, Xb, type="probs")[, which(levels(ytr)==pos)]
    } else if (m == "xgb") {
      d <- xgb.DMatrix(Xtr, label=as.numeric(ytr==pos))
      f <- xgb.train(list(objective="binary:logistic", nthread=1), d, nrounds=100, verbose=0)
      predict(f, xgb.DMatrix(Xte))
    } else if (m == "logistic") {

      cvfit <- cv.glmnet(Xtr_s, ytr, family="binomial", alpha=0, nfolds=5)
      as.numeric(predict(cvfit, Xte_s, s="lambda.min", type="response"))
    }
  }
  cv_roc <- function(X, y, groups, model, pos, grid) {
    ug <- unique(groups); tpr <- matrix(NA, n_iter, length(grid))
    for (i in seq_len(n_iter)) {
      set.seed(i); tr_g <- sample(ug, floor(0.8*length(ug)))
      tr <- groups %in% tr_g; te <- !tr
      if (length(unique(y[te]))<2 || length(unique(y[tr]))<2) next
      pv <- apply(X[tr,,drop=FALSE], 2, function(v) tryCatch(wilcox.test(v~y[tr])$p.value, error=function(e)1))
      sel <- which(p.adjust(pv,"fdr")<0.05); if (length(sel)<5) sel <- order(pv)[1:min(30,length(pv))]
      pr <- tryCatch(fit_prob(model, X[tr,sel,drop=FALSE], y[tr], X[te,sel,drop=FALSE], pos),
                     error=function(e) NULL)
      if (is.null(pr)) next
      r <- pROC::roc(as.integer(y[te]==pos), pr, quiet=TRUE, direction="<")
      co <- pROC::coords(r, "all", ret=c("specificity","sensitivity"), transpose=FALSE)
      fpr_i <- 1 - co$specificity; tpr_i <- co$sensitivity; o <- order(fpr_i)
      tpr[i, ] <- approx(fpr_i[o], tpr_i[o], grid, ties=max, rule=2)$y
    }
    data.table(fpr=grid, mean=colMeans(tpr,na.rm=TRUE), sd=apply(tpr,2,sd,na.rm=TRUE),
               auc=NA_real_)
  }
  load_X <- function(which) {
    if (which=="peak-list") {

      X <- readRDS(file.path(TXT_OBJ,"txt_matrix_aggregated.rds"))
      meta <- as.data.table(readRDS(file.path(TXT_OBJ,"txt_metadata_aggregated.rds")))
    } else {
      ds <- readRDS(file.path(RAW_OBJ,"raw_peak_strain_treatment_dataset.rds"))
      X <- ds$matrix; meta <- as.data.table(ds$metadata)
    }
    list(X=X, meta=meta)
  }

  spec <- list(
    list(ds="peak-list", task="res_sen", models=c("svm","xgb","adaboost"), tag="A", lab="peak-list — CRKp/CSKp"),
    list(ds="peak-list", task="kpc_ndm", models=c("rf","extra_trees","logistic"), tag="B", lab="peak-list — KPC/NDM"),
    list(ds="RAW dataset", task="res_sen", models=c("svm","xgb","rf"), tag="C", lab="RAW dataset — CRKp/CSKp"),
    list(ds="RAW dataset", task="kpc_ndm", models=c("adaboost","rf","xgb"), tag="D", lab="RAW dataset — KPC/NDM")
  )
  grid <- seq(0, 1, length.out=101)
  MLAB <- c(svm="SVM", rf="RF", xgb="XGB", adaboost="AdaBoost", extra_trees="Extra Trees", logistic="Logistic")
  panels <- list()
  for (s in spec) {
    dd <- load_X(s$ds); X <- dd$X; meta <- dd$meta
    if (s$task=="res_sen") { y <- factor(meta$label_res_sen, c("SEN","RES")); pos <- "RES" }
    else { keep <- meta$mechanism %in% c("KPC","NDM"); X <- X[keep,,drop=FALSE]
           y <- factor(meta$mechanism[keep], c("NDM","KPC")); pos <- "KPC"; meta <- meta[keep] }
    grp <- meta$strain_id
    res <- rbindlist(lapply(s$models, function(m) {
      cat(sprintf("  [%s/%s] %s\n", s$ds, s$task, m))
      r <- cv_roc(X, y, grp, m, pos, grid)
      r <- rbind(data.table(fpr=0, mean=0, sd=0, auc=NA_real_), r)
      r[, model := MLAB[m]][]
    }))
    panels[[s$tag]] <- ggplot(res, aes(fpr, mean, colour=model, fill=model)) +
      geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey60") +
      geom_ribbon(aes(ymin=pmax(0,mean-sd), ymax=pmin(1,mean+sd)), alpha=0.15, colour=NA) +
      geom_line(linewidth=0.7) +
      scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
      coord_fixed(ratio=1, xlim=c(0,1), ylim=c(0,1)) +
      labs(x="False positive rate", y="True positive rate", tag=s$tag, subtitle=s$lab, colour=NULL, fill=NULL) +
      theme_paper(11) + theme(legend.position=c(0.98,0.02), legend.justification=c(1,0),
                              plot.subtitle=element_text(size=10,face="bold"),
                              legend.background=element_rect(fill=alpha("white",0.7), colour=NA))
  }
  p <- (panels$A | panels$B) / (panels$C | panels$D)
  save_fig(p, "sfig5_roc_curves_suplementario", width = 24, height = 24)
}

if (do_fig(1)) fig_s1()
if (do_fig(2)) fig_s2()
if (do_fig(5)) fig_s5()
if (do_fig(3)) fig_s3()
if (do_fig(4)) fig_s4()
if (do_fig(7)) fig_s7()
if (do_fig(6)) fig_s6()
if (do_fig(8)) fig_s8()
cat("Done.\n")
