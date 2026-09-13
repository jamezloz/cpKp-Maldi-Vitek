source("scripts_txt/00_shared_utils.R")

suppressPackageStartupMessages({
  library(pROC); library(data.table)
  library(e1071); library(class); library(rpart); library(randomForest)
  library(ada); library(xgboost); library(glmnet)
})

sample_group_split <- function(groups, train_frac = 0.8) {
  ug <- unique(groups)
  tr_g <- sample(ug, size = floor(train_frac * length(ug)))
  list(train = groups %in% tr_g, test = !groups %in% tr_g)
}

fit_one_model <- function(model_name, Xtr, Xtr_sc, Xte, Xte_sc, ytr, yte, pos_label, neg_label, y_levels) {
  probs <- NULL; pred <- NULL
  if (model_name == "svm") {
    fit  <- svm(Xtr_sc, ytr, kernel = "radial", probability = TRUE)
    pp   <- predict(fit, Xte_sc, probability = TRUE)
    pred <- pp; probs <- attr(pp, "probabilities")[, pos_label]
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
    pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = y_levels)
  } else if (model_name == "extra_trees") {
    fit   <- randomForest(x = Xtr, y = ytr, ntree = 500, mtry = ncol(Xtr), replace = FALSE)
    probs <- predict(fit, Xte, type = "prob")[, pos_label]
    pred  <- predict(fit, Xte, type = "response")
  } else if (model_name == "adaboost") {
    Xtr_ada <- as.data.frame(Xtr); colnames(Xtr_ada) <- paste0("f", colnames(Xtr_ada))
    Xte_ada <- as.data.frame(Xte); colnames(Xte_ada) <- paste0("f", colnames(Xte_ada))
    fit       <- ada(x = Xtr_ada, y = ytr, iter = 100, loss = "exponential")
    probs_mat <- predict(fit, Xte_ada, type = "probs")
    probs     <- probs_mat[, which(levels(ytr) == pos_label)]
    pred      <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = y_levels)
  } else if (model_name == "logistic") {
    cvfit <- cv.glmnet(Xtr_sc, ytr, family = "binomial", alpha = 0, nfolds = 5)
    probs <- as.numeric(predict(cvfit, Xte_sc, s = "lambda.min", type = "response"))
    pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = y_levels)
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
    fit  <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", nthread = 1),
                      data = dtr, nrounds = 100, verbose = 0)
    probs <- predict(fit, dte)
    pred  <- factor(ifelse(probs >= 0.5, pos_label, neg_label), levels = y_levels)
  }
  list(probs = probs, pred = pred)
}

cv_global_featsel <- function(X_global, y, groups, model_name, n_iter = 100) {
  pos_label <- levels(y)[2]; neg_label <- levels(y)[1]; y_levels <- levels(y)
  rows <- list()
  for (i in seq_len(n_iter)) {
    sp  <- sample_group_split(groups)
    tr  <- sp$train; te <- sp$test
    ytr <- y[tr]; yte <- y[te]
    if (length(unique(ytr)) < 2 || length(unique(yte)) < 2) next

    Xtr <- X_global[tr, , drop = FALSE]
    Xte <- X_global[te, , drop = FALSE]
    Xtr_sc <- scale(Xtr)
    sc_ctr <- attr(Xtr_sc, "scaled:center"); sc_scl <- attr(Xtr_sc, "scaled:scale"); sc_scl[sc_scl == 0] <- 1
    Xte_sc <- scale(Xte, center = sc_ctr, scale = sc_scl)

    r <- tryCatch(fit_one_model(model_name, Xtr, Xtr_sc, Xte, Xte_sc, ytr, yte, pos_label, neg_label, y_levels),
                  error = function(e) { message(sprintf("  [iter %d %s error]: %s", i, model_name, e$message)); NULL })
    if (is.null(r) || is.null(r$pred) || is.null(r$probs)) next
    pred <- r$pred; probs <- r$probs

    TP <- sum(pred == pos_label & yte == pos_label); FN <- sum(pred == neg_label & yte == pos_label)
    TN <- sum(pred == neg_label & yte == neg_label); FP <- sum(pred == pos_label & yte == neg_label)
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    ba   <- mean(c(sens, spec), na.rm = TRUE)
    auc  <- tryCatch(as.numeric(pROC::auc(pROC::roc(as.numeric(yte == pos_label), probs, quiet = TRUE))),
                     error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.table(iter = i, model = model_name, accuracy = mean(pred == yte),
                                           balanced_accuracy = ba, auc = auc, sensitivity = sens,
                                           specificity = spec, PPV = ppv, n_features = ncol(X_global))
  }
  rbindlist(rows)
}

main <- function() {
  meta <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_metadata_aggregated.rds"))
  reps <- readRDS(file.path(TXT_MODULAR_OBJ, "txt_representations.rds"))
  X    <- reps$intensity_raw

  models <- c("rf", "xgb", "svm", "knn", "dt", "extra_trees", "adaboost", "logistic", "naive_bayes")

  tasks <- list(
    res_sen = list(keep = meta$label_res_sen %in% c("RES", "SEN"),
                   y_fn = function(m) factor(m$label_res_sen, levels = c("SEN", "RES"))),
    kpc_ndm = list(keep = meta$mechanism %in% c("KPC", "NDM"),
                   y_fn = function(m) factor(m$mechanism, levels = c("NDM", "KPC")))
  )

  all_res <- list()
  n_feat_used <- list()

  for (task in names(tasks)) {
    keep <- tasks[[task]]$keep
    Xk   <- X[keep, , drop = FALSE]
    y    <- tasks[[task]]$y_fn(meta[keep])
    grps <- meta$strain_id[keep]

    pvals <- apply(Xk, 2, function(v) tryCatch(wilcox.test(v ~ y)$p.value, error = function(e) 1))
    fdr   <- p.adjust(pvals, method = "fdr")
    keep_feat <- which(fdr < 0.05)
    if (length(keep_feat) == 0) keep_feat <- order(pvals)[seq_len(min(10, ncol(Xk)))]
    Xk_sel <- Xk[, keep_feat, drop = FALSE]
    n_feat_used[[task]] <- ncol(Xk_sel)
    cat(sprintf("Task %s: n=%d, features preseleccionadas globalmente=%d\n", task, sum(keep), ncol(Xk_sel)))

    for (mod in models) {
      cat("  model:", mod, "\n")
      res <- cv_global_featsel(Xk_sel, y, grps, model_name = mod, n_iter = 100)
      if (nrow(res)) { res[, task := task]; all_res[[paste(task, mod)]] <- res }
    }
  }

  out <- rbindlist(all_res, fill = TRUE)
  save_table(out, "ml_global_featsel_iterations")

  summary_dt <- out[, .(
    accuracy_mean = mean(accuracy, na.rm = TRUE), balanced_acc_mean = mean(balanced_accuracy, na.rm = TRUE),
    balanced_acc_sd = sd(balanced_accuracy, na.rm = TRUE), auc_mean = mean(auc, na.rm = TRUE), auc_sd = sd(auc, na.rm = TRUE),
    sensitivity_mean = mean(sensitivity, na.rm = TRUE), specificity_mean = mean(specificity, na.rm = TRUE),
    PPV_mean = mean(PPV, na.rm = TRUE), n_features = n_feat_used[[task[1]]], n_iter = sum(!is.na(balanced_accuracy))
  ), by = .(task, model)]
  save_table(summary_dt, "ml_global_featsel_summary")

  cat("\nResumen (preseleccion global, 9 modelos):\n")
  print(summary_dt[order(task, -balanced_acc_mean)])
  cat("Saved to:", TXT_MODULAR_TAB, "\n")
}

main()
