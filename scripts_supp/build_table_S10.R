suppressPackageStartupMessages(library(data.table))
OUT <- "resultados/PAPER_SUPP_TABLES"
TXT_TAB <- "resultados/TXT_modular/tables"
RAW_TAB <- "resultados/RAW_modular/tables"

lab <- c(rf="RF", xgb="XGB", svm="SVM", knn="KNN", dt="DT",
        extra_trees="ET", adaboost="ADA", logistic="LOG", naive_bayes="NB")
task_lab <- c(res_sen="CRKp/CSKp", kpc_ndm="KPC/NDM")

intra <- fread(file.path(TXT_TAB, "ml_extended_classifiers_summary.csv"))
intra <- intra[, .(task = task_lab[task], model = lab[model], strategy = "Intra-fold Wilcoxon",
                   BA = balanced_acc_mean, BA_sd = balanced_acc_sd, n_features = "within-fold")]

glob <- fread(file.path(TXT_TAB, "ml_global_featsel_summary.csv"))
glob <- glob[, .(task = task_lab[task], model = lab[model], strategy = "Global preselection",
                 BA = balanced_acc_mean, BA_sd = balanced_acc_sd, n_features = as.character(n_features))]

nofs <- fread(file.path(TXT_TAB, "ml_extended_nofeatsel_summary.csv"))
nofs <- nofs[, .(task = task_lab[task], model = lab[model], strategy = "No selection (all features)",
                 BA = balanced_acc_mean, BA_sd = balanced_acc_sd, n_features = "3,958")]

peaklist <- rbindlist(list(intra, glob, nofs))
fwrite(peaklist, file.path(OUT, "S10_peaklist_featsel_comparison.csv"))

raw <- fread(file.path(RAW_TAB, "raw_extended_classifiers_strain_treatment.csv"))
raw_sum <- raw[, .(BA = mean(balanced_accuracy), BA_sd = sd(balanced_accuracy)), by = .(task, model)]
raw_sum <- raw_sum[, .(task = task_lab[task], model = lab[model], strategy = "Intra-fold Wilcoxon (only strategy run)",
                       BA = BA, BA_sd = BA_sd, n_features = "within-fold")]
fwrite(raw_sum, file.path(OUT, "S10_raw_intrafold.csv"))

cat("Peak-list rows:", nrow(peaklist), "| RAW rows:", nrow(raw_sum), "\n")
cat("\n--- feature counts used ---\n")
print(unique(glob[, .(task, n_features)]))
