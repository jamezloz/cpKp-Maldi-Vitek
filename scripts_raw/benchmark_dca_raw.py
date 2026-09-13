from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupShuffleSplit

RANDOM_STATE = 42
N_SPLITS = 30
PT_GRID = np.linspace(0, 0.99, 200)

TABLE_DIR = Path("resultados/RAW_modular/tables")
OUT_DIR = TABLE_DIR

def compute_net_benefit(y_true: np.ndarray, y_prob: np.ndarray, pt_grid: np.ndarray) -> np.ndarray:
    n = len(y_true)
    nb = np.empty(len(pt_grid))
    for i, pt in enumerate(pt_grid):
        pred = (y_prob >= pt).astype(int)
        tp = int(((pred == 1) & (y_true == 1)).sum())
        fp = int(((pred == 1) & (y_true == 0)).sum())
        if pt >= 1.0:
            nb[i] = 0.0
        else:
            nb[i] = tp / n - fp / n * (pt / (1.0 - pt))
    return nb

def treat_all_nb(y_true: np.ndarray, pt_grid: np.ndarray) -> np.ndarray:
    n = len(y_true)
    n_pos = y_true.sum()
    n_neg = n - n_pos
    nb = np.empty(len(pt_grid))
    for i, pt in enumerate(pt_grid):
        if pt >= 1.0:
            nb[i] = 0.0
        else:
            nb[i] = n_pos / n - n_neg / n * (pt / (1.0 - pt))
    return nb

def run_dca(task: str, Xmat: np.ndarray, y: np.ndarray, strain_ids: np.ndarray) -> pd.DataFrame:
    print(f"\n=== {task} ===")
    gss = GroupShuffleSplit(n_splits=N_SPLITS, test_size=0.1, random_state=RANDOM_STATE)
    all_y_true, all_y_prob = [], []

    for train_idx, test_idx in gss.split(Xmat, y, groups=strain_ids):
        Xtr, Xte = Xmat[train_idx], Xmat[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        if len(np.unique(ytr)) < 2:
            continue
        rf = RandomForestClassifier(
            n_estimators=200,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        rf.fit(Xtr, ytr)
        prob = rf.predict_proba(Xte)[:, 1]
        all_y_true.extend(yte.tolist())
        all_y_prob.extend(prob.tolist())

    y_true_all = np.array(all_y_true)
    y_prob_all = np.array(all_y_prob)
    print(f"  Collected {len(y_true_all)} predictions. Prevalence: {y_true_all.mean():.3f}")

    nb_rf = compute_net_benefit(y_true_all, y_prob_all, PT_GRID)
    nb_treat_all = treat_all_nb(y_true_all, PT_GRID)

    return pd.DataFrame({
        "threshold": PT_GRID,
        "net_benefit_rf": nb_rf,
        "net_benefit_treat_all": nb_treat_all,
        "net_benefit_treat_none": np.zeros(len(PT_GRID)),
    })

print("Loading strain×treatment dataset (corrected, 539 obs)...")
X_df = pd.read_csv(TABLE_DIR / "raw_peak_strain_treatment_matrix.csv", index_col=0)
meta = pd.read_csv(TABLE_DIR / "raw_peak_strain_treatment_metadata.csv")
meta.index = meta["sample"]

X_df["_label_res_sen"] = meta.loc[X_df.index, "label_res_sen"].values
X_df["_label_mech"]    = meta.loc[X_df.index, "mechanism"].values
X_df["_strain_raw"]    = meta.loc[X_df.index, "strain_id"].values

X_df = X_df.dropna(subset=["_label_res_sen"])
print(f"  After label join: {X_df.shape[0]} rows")

feature_cols_all = [c for c in X_df.columns if not c.startswith("_")]
kruskal = pd.read_csv(TABLE_DIR / "raw_classic_biomarkers_kruskal_global.csv")
canon = set(kruskal["feature_chr"].astype(str))
feature_cols = [c for c in feature_cols_all if str(c) in canon]
if len(feature_cols) != len(canon):
    raise SystemExit(f"Features de R no presentes en la matriz: "
                     f"{len(canon)} esperadas, {len(feature_cols)} encontradas")
print(f"  Features (prevalencia maxima por clase >= 5%, via R): {len(feature_cols)}")
Xmat = X_df[feature_cols].values.astype(float)
strain_ids = X_df["_strain_raw"].values

y_res_sen = (X_df["_label_res_sen"] == "RES").astype(int).values
df_dca_res = run_dca("res_sen", Xmat, y_res_sen, strain_ids)
out_res = OUT_DIR / "dca_res_sen_raw.csv"
df_dca_res.to_csv(out_res, index=False)
print(f"  Saved: {out_res}")

mask_kpc_ndm = X_df["_label_mech"].isin(["KPC", "NDM"])
X_kpc = Xmat[mask_kpc_ndm]
y_kpc = (X_df.loc[mask_kpc_ndm, "_label_mech"] == "KPC").astype(int).values
s_kpc = strain_ids[mask_kpc_ndm]
print(f"\n  kpc_ndm subset: {X_kpc.shape[0]} rows, KPC={y_kpc.sum()}, NDM={(y_kpc==0).sum()}")
df_dca_kpc = run_dca("kpc_ndm", X_kpc, y_kpc, s_kpc)
out_kpc = OUT_DIR / "dca_kpc_ndm_raw.csv"
df_dca_kpc.to_csv(out_kpc, index=False)
print(f"  Saved: {out_kpc}")

print("\nPaper 1 DCA complete.")
