from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupShuffleSplit

RANDOM_STATE = 42
N_SPLITS     = 30
PT_GRID      = np.linspace(0, 0.99, 200)

TABLE_DIR = Path("resultados/TXT_modular/tables")
OUT_DIR   = TABLE_DIR

def compute_net_benefit(y_true: np.ndarray, y_prob: np.ndarray,
                        pt_grid: np.ndarray) -> np.ndarray:
    n  = len(y_true)
    nb = np.empty(len(pt_grid))
    for i, pt in enumerate(pt_grid):
        pred = (y_prob >= pt).astype(int)
        tp   = int(((pred == 1) & (y_true == 1)).sum())
        fp   = int(((pred == 1) & (y_true == 0)).sum())
        nb[i] = 0.0 if pt >= 1.0 else tp / n - fp / n * (pt / (1.0 - pt))
    return nb

def treat_all_nb(y_true: np.ndarray, pt_grid: np.ndarray) -> np.ndarray:
    n, n_pos, n_neg = len(y_true), y_true.sum(), (y_true == 0).sum()
    nb = np.empty(len(pt_grid))
    for i, pt in enumerate(pt_grid):
        nb[i] = 0.0 if pt >= 1.0 else n_pos / n - n_neg / n * (pt / (1.0 - pt))
    return nb

def run_dca(task: str, Xmat: np.ndarray, y: np.ndarray,
            strain_ids: np.ndarray) -> pd.DataFrame:
    print(f"\n=== {task} ===")
    gss = GroupShuffleSplit(n_splits=N_SPLITS, test_size=0.1, random_state=RANDOM_STATE)
    all_y_true, all_y_prob = [], []

    for train_idx, test_idx in gss.split(Xmat, y, groups=strain_ids):
        Xtr, Xte = Xmat[train_idx], Xmat[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        if len(np.unique(ytr)) < 2:
            continue
        rf = RandomForestClassifier(n_estimators=200, class_weight="balanced",
                                    random_state=RANDOM_STATE, n_jobs=-1)
        rf.fit(Xtr, ytr)
        prob = rf.predict_proba(Xte)[:, 1]
        all_y_true.extend(yte.tolist())
        all_y_prob.extend(prob.tolist())

    y_true_all = np.array(all_y_true)
    y_prob_all = np.array(all_y_prob)
    print(f"  Collected {len(y_true_all)} predictions. Prevalence: {y_true_all.mean():.3f}")

    nb_rf        = compute_net_benefit(y_true_all, y_prob_all, PT_GRID)
    nb_treat_all = treat_all_nb(y_true_all, PT_GRID)

    return pd.DataFrame({
        "threshold":              PT_GRID,
        "net_benefit_rf":         nb_rf,
        "net_benefit_treat_all":  nb_treat_all,
        "net_benefit_treat_none": np.zeros(len(PT_GRID)),
    })

print("Loading TXT matrix and metadata...")
X_df = pd.read_csv(TABLE_DIR / "txt_matrix_aggregated.csv", index_col=0)
meta = pd.read_csv(TABLE_DIR / "txt_metadata_aggregated.csv")

kruskal = pd.read_csv(TABLE_DIR / "biomarkers_global_kruskal.csv")
canon = set(kruskal["feature"].astype(str))
feat_cols_f = [c for c in X_df.columns if str(c) in canon]
if len(feat_cols_f) != len(canon):
    raise SystemExit(f"Features de R no presentes en la matriz: "
                     f"{len(canon)} esperadas, {len(feat_cols_f)} encontradas")
print(f"  Features (prevalencia maxima por clase >= 5%, via R): {len(feat_cols_f)}")

Xmat       = X_df[feat_cols_f].values.astype(float)
meta["strain_id"] = meta["sample"].str.split("__").str[0]
strain_ids = meta["strain_id"].values

mask_rs  = meta["label_res_sen"].isin(["RES", "SEN"])
y_res    = (meta.loc[mask_rs, "label_res_sen"] == "RES").astype(int).values
s_res    = strain_ids[mask_rs.values]
df_res   = run_dca("res_sen", Xmat[mask_rs.values], y_res, s_res)
out_res  = OUT_DIR / "dca_res_sen_txt.csv"
df_res.to_csv(out_res, index=False)
print(f"  Saved: {out_res}")

mask_kn  = meta["mechanism"].isin(["KPC", "NDM"])
y_kn     = (meta.loc[mask_kn, "mechanism"] == "KPC").astype(int).values
s_kn     = strain_ids[mask_kn.values]
df_kn    = run_dca("kpc_ndm", Xmat[mask_kn.values], y_kn, s_kn)
out_kn   = OUT_DIR / "dca_kpc_ndm_txt.csv"
df_kn.to_csv(out_kn, index=False)
print(f"  Saved: {out_kn}")

print("\nPaper 1 TXT DCA complete.")
