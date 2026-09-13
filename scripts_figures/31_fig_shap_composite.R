source("scripts_figures/00_fig_utils.R")

shap_rs  <- fread(file.path(TXT_TAB, "shap_matrix_res_sen.csv"))
imp_rs   <- fread(file.path(TXT_TAB, "shap_importance_res_sen.csv"))
shap_kpc <- fread(file.path(TXT_TAB, "shap_matrix_kpc_ndm.csv"))
imp_kpc  <- fread(file.path(TXT_TAB, "shap_importance_kpc_ndm.csv"))
stab_rs  <- fread(file.path(TXT_TAB, "stability_selection_res_sen.csv"))
stab_kpc <- fread(file.path(TXT_TAB, "stability_selection_kpc_ndm.csv"))
raw_imp_rs  <- fread(file.path(RAW_TAB, "shap_importance_res_sen_raw.csv"))
raw_imp_kpc <- fread(file.path(RAW_TAB, "shap_importance_kpc_ndm_raw.csv"))

build_beeswarm <- function(shap_mat, imp_dt, n_top = 20) {
  top_feats <- as.character(imp_dt$feature[1:min(n_top, nrow(imp_dt))])
  meta_cols  <- c("true_label", "pred_prob")
  feat_cols  <- setdiff(names(shap_mat), meta_cols)
  feat_cols  <- feat_cols[feat_cols %in% top_feats]

  shap_long <- melt(shap_mat[, c(feat_cols, "true_label"), with = FALSE],
                    id.vars = "true_label",
                    variable.name = "feature",
                    value.name = "shap_value")
  shap_long[, feature := as.integer(as.character(feature))]
  shap_long <- merge(shap_long,
                     imp_dt[, .(feature, mean_abs_shap, direction)],
                     by = "feature", all.x = TRUE)
  set.seed(42)
  shap_long[, feat_label := paste0(feature, " Da")]
  shap_long[, feat_label := factor(feat_label,
    levels = paste0(rev(imp_dt$feature[1:n_top]), " Da"))]
  shap_long
}

stab_rs_feats  <- stab_rs$feature[stab_rs$selection_freq >= 1.0]
stab_kpc_feats <- stab_kpc$feature[stab_kpc$selection_freq >= 1.0]

bee_rs <- build_beeswarm(shap_rs, imp_rs, n_top = 20)
imp_rs[, stable := feature %in% stab_rs_feats]

pA <- ggplot(bee_rs, aes(shap_value, feat_label, colour = shap_value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_jitter(height = 0.28, size = 1.2, alpha = 0.60, width = 0) +
  scale_colour_gradient2(low = "#2980B9", mid = "#F0F0F0", high = "#C0392B",
                         midpoint = 0,
                         name = "SHAP\n(blue=CSKp, red=CRKp)") +
  labs(x = "SHAP value", y = NULL, tag = "A",
       subtitle = "TXT · CRKp vs CSKp · top 20") +
  theme_paper() +
  theme(legend.position  = "right",
        axis.text.y      = element_text(size = 8.5),
        plot.subtitle    = element_text(size = 10, face = "bold"))

bee_kpc <- build_beeswarm(shap_kpc, imp_kpc, n_top = 15)

pB <- ggplot(bee_kpc, aes(shap_value, feat_label, colour = shap_value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_jitter(height = 0.28, size = 1.2, alpha = 0.60, width = 0) +
  scale_colour_gradient2(low = "#27AE60", mid = "#F0F0F0", high = "#E67E22",
                         midpoint = 0,
                         name = "SHAP\n(green=NDM, orange=KPC)") +
  labs(x = "SHAP value", y = NULL, tag = "B",
       subtitle = "TXT · KPC vs NDM · top 15") +
  theme_paper() +
  theme(legend.position  = "right",
        axis.text.y      = element_text(size = 8.5),
        plot.subtitle    = element_text(size = 10, face = "bold"))

raw_imp_rs_top <- raw_imp_rs[1:min(15, .N)]
raw_imp_rs_top[, feat_label := paste0(round(feature, 0), " Da")]
raw_imp_rs_top[, feat_label := factor(feat_label, levels = rev(feat_label))]

pC <- ggplot(raw_imp_rs_top, aes(mean_abs_shap, feat_label)) +
  geom_col(fill = "#C0392B", alpha = 0.80, width = 0.72) +
  labs(x = "Mean |SHAP|", y = NULL, tag = "C",
       subtitle = "RAW · CRKp vs CSKp · top 15") +
  theme_paper() +
  theme(axis.text.y   = element_text(size = 8.5),
        plot.subtitle = element_text(size = 10, face = "bold"))

raw_imp_kpc_top <- raw_imp_kpc[1:min(15, .N)]
raw_imp_kpc_top[, feat_label := paste0(round(feature, 0), " Da")]
raw_imp_kpc_top[, feat_label := factor(feat_label, levels = rev(feat_label))]

pD <- ggplot(raw_imp_kpc_top, aes(mean_abs_shap, feat_label)) +
  geom_col(fill = "#E67E22", alpha = 0.80, width = 0.72) +
  labs(x = "Mean |SHAP|", y = NULL, tag = "D",
       subtitle = "RAW · KPC vs NDM · top 15") +
  theme_paper() +
  theme(axis.text.y   = element_text(size = 8.5),
        plot.subtitle = element_text(size = 10, face = "bold"))

fig_out <- (pA | pB) / (pC | pD)
save_fig(fig_out, "fig_shap_composite", width = 34, height = 32)
cat("Saved: fig_shap_composite.png\n")
