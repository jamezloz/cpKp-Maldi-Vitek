suppressPackageStartupMessages(library(data.table))
OUT <- "resultados/PAPER_SUPP_TABLES"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TXT_TAB <- "resultados/TXT_modular/tables"
RAW_TAB <- "resultados/RAW_modular/tables"

lab <- c(rf="RF", xgb="XGB", svm="SVM", knn="KNN", dt="DT",
        extra_trees="ET", adaboost="ADA", logistic="LOG", naive_bayes="NB")
task_lab <- c(res_sen="CRKp/CSKp", kpc_ndm="KPC/NDM")

a1 <- fread(file.path(TXT_TAB, "ml_extended_classifiers_summary.csv"))
a1 <- a1[, .(task = task_lab[task], level = "strain x protocol", model = lab[model],
            BA = balanced_acc_mean, BA_sd = balanced_acc_sd,
            AUC = auc_mean, AUC_sd = auc_sd,
            sens = sensitivity_mean, spec = specificity_mean, ppv = PPV_mean)]

rl <- fread(file.path(TXT_TAB, "ml_raw_level_iterations.csv"))
a2 <- rl[, .(BA = mean(balanced_accuracy), BA_sd = sd(balanced_accuracy),
            AUC = mean(auc), AUC_sd = sd(auc)), by = .(task, model)]
a2[, task := task_lab[task]]; a2[, model := lab[model]]
a2[, level := "replicate"]; a2[, sens := NA_real_]; a2[, spec := NA_real_]; a2[, ppv := NA_real_]
setcolorder(a2, names(a1))

partA <- rbindlist(list(a1, a2))
setorder(partA, task, level, -BA)
fwrite(partA, file.path(OUT, "S8_partA_peaklist.csv"))
cat("Part A (peak-list): strain x protocol =", nrow(a1), "rows | replicate =", nrow(a2), "rows\n")

b1 <- fread(file.path(RAW_TAB, "raw_extended_classifiers_strain_treatment.csv"))
b1s <- b1[, .(BA = mean(balanced_accuracy), BA_sd = sd(balanced_accuracy),
             AUC = mean(auc), AUC_sd = sd(auc),
             sens = mean(sensitivity), spec = mean(specificity)), by = .(task, model)]
b1s[, task := task_lab[task]]; b1s[, model := lab[model]]; b1s[, level := "strain x protocol"]

b2 <- fread(file.path(RAW_TAB, "levels_comparison/ml_levels_iterations_replicate100.csv"))
b2s <- b2[, .(BA = mean(balanced_accuracy), BA_sd = sd(balanced_accuracy),
             AUC = mean(auc), AUC_sd = sd(auc),
             sens = mean(sensitivity), spec = mean(specificity)), by = .(task, model)]
b2s[, task := task_lab[task]]; b2s[, model := lab[model]]; b2s[, level := "replicate"]

partB <- rbindlist(list(b1s, b2s), fill = TRUE)
partB[, ppv := NA_real_]
setcolorder(partB, names(partA))
setorder(partB, task, level, -BA)
fwrite(partB, file.path(OUT, "S8_partB_raw.csv"))
cat("Part B (RAW): strain x protocol =", nrow(b1s), "rows | replicate =", nrow(b2s), "rows\n")

cat("\n--- Sanity: peak-list N (replicate) ---\n")
meta <- readRDS("resultados/TXT_modular/objects/txt_metadata_raw.rds")
cat("TXT replicate total spectra:", nrow(meta), "\n")
