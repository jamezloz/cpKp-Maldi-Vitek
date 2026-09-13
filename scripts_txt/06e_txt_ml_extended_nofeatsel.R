source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(pROC)
  library(data.table)
  library(e1071)
  library(class)
  library(rpart)
  library(randomForest)
  library(ada)
  library(xgboost)
  library(glmnet)
})

sample_group_split <- function(groups, train_frac = 0.8) {
  ug <- unique(groups)
  tr_g <- sample(ug, size = floor(train_frac * length(ug)))
  list(train = groups %in% tr_g, test = !groups %in% tr_g)
}

cv_extended_nofeatsel <- function(X, y, groups, model_name, n_iter = 100) {
  pos_label <- levels(y)[2]
  neg_label <- levels(y)[1]
  rows <- list()

  for (i in seq_len(n_iter)) {
    sp  <- sample_group_split(groups)
    tr  <- sp$train; te  <- sp$test
    ytr <- y[tr];    yte <- y[te]
    if (length(unique(ytr)) < 2 || length(unique(yte)) < 2) next

    Xtr_full <- X[tr, , drop = FALSE]
    Xte_full <- X[te, , drop = FALSE]
    col_var  <- apply(Xtr_full, 2, var)
    keep_var <- which(col_var > 0 & is.finite(col_var))
    Xtr <- Xtr_full[, keep_var, drop = FALSE]
    Xte <- Xte_full[, keep_var, drop = FALSE]

    Xtr_sc <- scale(Xtr)
    sc_ctr  <- attr(Xtr_sc, "scaled:center")
    sc_scl  <- attr(Xtr_sc, "scaled:scale")
    sc_scl[sc_scl == 0] <- 1
    Xte_sc  <- scale(Xte, center = sc_ctr, scale = sc_scl)

    probs <- NULL
    pred  <- NULL

    tryCatch({
      if (model_name == "svm") {
        fit   <- svm(Xtr_sc, ytr, kernel = "radial", probability = TRUE)
        pp    <- predict(fit, Xte_sc, probability = TRUE)
        pred  <- pp
        probs <- attr(pp, "probabilities")[, pos_label]

      } else if (model_name == "knn") {
        k     <- max(1, round(sqrt(nrow(Xtr_sc))))
        pred  <- knn(Xtr_sc, Xte_sc, cl = ytr, k = k, prob = TRUE)
        raw_p <- attr(pred, "prob")
        probs <- ifelse(pred == pos_label, raw_p, 1 - raw_p)

      } else if (model_name == "dt") {
        fit   <- rpart(y ~ ., data = data.frame(y = ytr, Xtr), method = "class",
                       control = rpart.control(minsplit = 5, cp = 0.001))
        pp    <- predict(fit, data.frame(Xte), type = "prob")
        probs <- pp[, pos_label]
        pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = levels(y))

      } else if (model_name == "extra_trees") {
        fit   <- randomForest(x = Xtr, y = ytr, ntree = 500,
                              mtry = ncol(Xtr), replace = FALSE)
        probs <- predict(fit, Xte, type = "prob")[, pos_label]
        pred  <- predict(fit, Xte, type = "response")

      } else if (model_name == "adaboost") {
        Xtr_ada <- as.data.frame(Xtr); colnames(Xtr_ada) <- paste0("f", colnames(Xtr_ada))
        Xte_ada <- as.data.frame(Xte); colnames(Xte_ada) <- paste0("f", colnames(Xte_ada))
        fit        <- ada(x = Xtr_ada, y = ytr, iter = 100, loss = "exponential")
        probs_mat  <- predict(fit, Xte_ada, type = "probs")
        probs      <- probs_mat[, which(levels(ytr) == pos_label)]
        pred       <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = levels(y))

      } else if (model_name == "logistic") {

        cvfit <- cv.glmnet(Xtr_sc, ytr, family = "binomial", alpha = 0, nfolds = 5)
        probs <- as.numeric(predict(cvfit, Xte_sc, s = "lambda.min", type = "response"))
        pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = levels(y))

      } else if (model_name == "naive_bayes") {
        fit   <- naiveBayes(Xtr_sc, ytr)
        pp    <- predict(fit, Xte_sc, type = "raw")
        probs <- pp[, pos_label]
        pred  <- predict(fit, Xte_sc, type = "class")

      } else if (model_name == "rf") {
        fit   <- randomForest(x = Xtr, y = ytr, ntree = 500)
        probs <- predict(fit, Xte, type = "prob")[, pos_label]
        pred  <- predict(fit, Xte, type = "response")

      } else if (model_name == "xgb") {
        dtr  <- xgb.DMatrix(data = Xtr, label = as.numeric(ytr == pos_label))
        dte  <- xgb.DMatrix(data = Xte, label = as.numeric(yte == pos_label))
        fit  <- xgb.train(
          params  = list(objective = "binary:logistic", eval_metric = "auc", nthread = 1),
          data    = dtr, nrounds = 100, verbose = 0
        )
        probs <- predict(fit, dte)
        pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = levels(y))
      }
    }, error = function(e) {
      message(sprintf("  [iter %d %s error]: %s", i, model_name, e$message))
    })

    if (is.null(pred) || is.null(probs)) next

    TP <- sum(pred == pos_label & yte == pos_label)
    FN <- sum(pred == neg_label & yte == pos_label)
    TN <- sum(pred == neg_label & yte == neg_label)
    FP <- sum(pred == pos_label & yte == neg_label)

    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    ba   <- mean(c(sens, spec), na.rm = TRUE)
    auc  <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(as.numeric(yte == pos_label), probs, quiet = TRUE))),
      error = function(e) NA_real_
    )

    rows[[length(rows) + 1]] <- data.table(
      iter = i, model = model_name,
      accuracy = mean(pred == yte),
      balanced_accuracy = ba, auc = auc,
      sensitivity = sens, specificity = spec, PPV = ppv,
      n_features = ncol(Xtr)
    )
  }
  rbindlist(rows)
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X    <- reps$intensity_raw
  grp  <- meta$strain_id

  bin      <- X > 0
  res_i    <- meta$label_res_sen == "RES"
  prev_max <- pmax(colMeans(bin[res_i, , drop = FALSE]),
                   colMeans(bin[!res_i, , drop = FALSE]))
  Xf       <- X[, prev_max >= 0.05, drop = FALSE]

  y_rs   <- factor(meta$label_res_sen, levels = c("SEN", "RES"))
  p_rs   <- apply(Xf, 2, function(v) tryCatch(wilcox.test(v ~ y_rs)$p.value, error = function(e) NA_real_))
  feats_rs <- colnames(Xf)[which(!is.na(p_rs) & p.adjust(p_rs, "fdr") < 0.05)]

  keep_kn <- meta$mechanism %in% c("KPC", "NDM")
  y_kn    <- factor(meta$mechanism[keep_kn], levels = c("NDM", "KPC"))
  p_kn    <- apply(Xf[keep_kn, , drop = FALSE], 2,
                   function(v) tryCatch(wilcox.test(v ~ y_kn)$p.value, error = function(e) NA_real_))
  feats_kn <- colnames(Xf)[which(!is.na(p_kn) & p.adjust(p_kn, "fdr") < 0.05)]
  cat(sprintf("RES/SEN global features: %d  |  KPC/NDM global features: %d\n",
              length(feats_rs), length(feats_kn)))

  models <- c("rf", "xgb", "svm", "knn", "dt", "extra_trees", "adaboost", "logistic", "naive_bayes")

  out_dir   <- file.path(TXT_MODULAR_OUT, "tables")
  iter_file <- file.path(out_dir, "ml_extended_nofeatsel_iterations.csv")
  if (file.exists(iter_file)) file.remove(iter_file)

  save_rows <- function(res, task_name) {
    if (nrow(res) == 0) return(invisible(NULL))
    res[, task := task_name]
    fwrite(res, iter_file, append = file.exists(iter_file))
    cat(sprintf("  -> %d rows saved\n", nrow(res)))
  }

  y_res  <- factor(meta$label_res_sen, levels = c("SEN", "RES"))
  X_rs   <- X[, feats_rs, drop = FALSE]
  for (mod in models) {
    cat(sprintf("[RES/SEN] %s ... ", mod)); flush.console()
    t0  <- proc.time()["elapsed"]
    res <- cv_extended_nofeatsel(X_rs, y_res, grp, model_name = mod, n_iter = 100)
    cat(sprintf("%.1fs\n", proc.time()["elapsed"] - t0))
    save_rows(res, "res_sen")
  }

  keep_mech <- meta$mechanism %in% c("KPC", "NDM")
  y_mech    <- factor(meta$mechanism[keep_mech], levels = c("NDM", "KPC"))
  X_kn      <- X[keep_mech, feats_kn, drop = FALSE]
  grp_m     <- grp[keep_mech]
  for (mod in models) {
    cat(sprintf("[KPC/NDM] %s ... ", mod)); flush.console()
    t0  <- proc.time()["elapsed"]
    res <- cv_extended_nofeatsel(X_kn, y_mech, grp_m, model_name = mod, n_iter = 100)
    cat(sprintf("%.1fs\n", proc.time()["elapsed"] - t0))
    save_rows(res, "kpc_ndm")
  }

  out <- fread(iter_file)

  summary_dt <- out[, .(
    accuracy_mean     = mean(accuracy,          na.rm = TRUE),
    balanced_acc_mean = mean(balanced_accuracy, na.rm = TRUE),
    balanced_acc_sd   = sd(balanced_accuracy,   na.rm = TRUE),
    auc_mean          = mean(auc,               na.rm = TRUE),
    auc_sd            = sd(auc,                 na.rm = TRUE),
    sensitivity_mean  = mean(sensitivity,       na.rm = TRUE),
    specificity_mean  = mean(specificity,       na.rm = TRUE),
    PPV_mean          = mean(PPV,               na.rm = TRUE),
    n_iter            = .N
  ), by = .(task, model)]

  fwrite(summary_dt, file.path(out_dir, "ml_extended_nofeatsel_summary.csv"))

  cat("\n=== SUMMARY (no feature selection) ===\n")
  print(summary_dt[order(task, -balanced_acc_mean)])
  cat("Saved to:", out_dir, "\n")
}

set.seed(42)
main()
