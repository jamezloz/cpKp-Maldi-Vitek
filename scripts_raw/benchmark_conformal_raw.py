from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import LabelEncoder

try:
    from mapie.classification import SplitConformalClassifier
except ImportError:
    raise SystemExit("pip install mapie")

RANDOM_STATE = 42
N_OUTER      = 30
ALPHAS       = [0.05, 0.10, 0.20]
TABLE_DIR    = Path("resultados/RAW_modular/tables")
OUT_DIR      = TABLE_DIR

def extract_strain(fname: str) -> str:
    import re
    base = fname.replace(".mzML", "").replace(".txt", "").replace(".TXT", "")
    m = re.split(r"[_\-](?:col|ext)[_\-]", base, maxsplit=1)
    if len(m) > 1:
        return m[0].lower()
    parts = base.rsplit("-", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0].lower()
    return base.lower()

def run_conformal(task: str, Xmat: np.ndarray, y: np.ndarray,
                  strain_ids: np.ndarray, n_outer: int = N_OUTER) -> pd.DataFrame:
    print(f"\n=== {task} ===")
    gss = GroupShuffleSplit(n_splits=n_outer, test_size=0.2, random_state=RANDOM_STATE)
    rows = []

    for fold_i, (train_idx, test_idx) in enumerate(gss.split(Xmat, y, groups=strain_ids)):
        X_tr, X_te = Xmat[train_idx], Xmat[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]

        strain_tr = strain_ids[train_idx]
        gss_cal = GroupShuffleSplit(n_splits=1, test_size=0.2,
                                    random_state=RANDOM_STATE + fold_i)
        prop_rel, cal_rel = next(gss_cal.split(X_tr, y_tr, groups=strain_tr))
        prop_idx = train_idx[prop_rel]
        cal_idx = train_idx[cal_rel]
        X_prop, X_cal = Xmat[prop_idx], Xmat[cal_idx]
        y_prop, y_cal = y[prop_idx], y[cal_idx]

        clf = RandomForestClassifier(n_estimators=500, random_state=RANDOM_STATE,
                                     n_jobs=-1, class_weight="balanced")
        clf.fit(X_prop, y_prop)

        results_conf = {}
        for alpha in ALPHAS:
            conf_level = 1 - alpha
            mapie = SplitConformalClassifier(
                clf, confidence_level=conf_level, prefit=True,
                random_state=RANDOM_STATE)
            mapie.conformalize(X_cal, y_cal)
            results_conf[conf_level] = mapie

        y_pred_point = clf.predict(X_te)
        ba = balanced_accuracy_score(y_te, y_pred_point)

        try:
            prob = clf.predict_proba(X_te)
            classes = list(clf.classes_)
            if len(classes) == 2:
                pos_i = classes.index(1) if 1 in classes else 1
                auc = roc_auc_score((y_te == (classes[1] if len(classes) > 1 else 1)).astype(int),
                                    prob[:, pos_i])
            else:
                auc = np.nan
        except Exception:
            auc = np.nan

        fold_row = {"fold": fold_i, "balanced_accuracy": ba, "auc": auc}

        le = LabelEncoder().fit(np.unique(y))
        y_te_enc = le.transform(y_te)
        for alpha in ALPHAS:
            conf_level = 1 - alpha
            mapie_c = results_conf[conf_level]

            _, sets = mapie_c.predict_set(X_te)
            sets_sq = sets[:, :, 0]
            set_sizes = sets_sq.sum(axis=1)
            covered = np.array([bool(sets_sq[i, y_te_enc[i]]) for i in range(len(y_te))])
            lvl = int(conf_level * 100)
            fold_row[f"coverage_{lvl}"] = float(covered.mean())
            fold_row[f"singleton_{lvl}"] = float((set_sizes == 1).mean())
            fold_row[f"empty_{lvl}"] = float((set_sizes == 0).mean())
            fold_row[f"mean_set_size_{lvl}"] = float(set_sizes.mean())
        rows.append(fold_row)

    df = pd.DataFrame(rows)
    print(f"  BA={df['balanced_accuracy'].mean():.4f}  AUC={df['auc'].mean():.4f}")
    for alpha in ALPHAS:
        lvl = int((1-alpha)*100)
        print(f"  α={alpha}: coverage={df[f'coverage_{lvl}'].mean():.4f}  "
              f"singleton={df[f'singleton_{lvl}'].mean():.4f}")
    return df.mean(numeric_only=True).to_frame().T.assign(task=task, model="rf", aggregation="replicate", n_folds=n_outer)

print("Loading strain×treatment dataset (corrected, 539 obs)...")
X_df = pd.read_csv(TABLE_DIR / "raw_peak_strain_treatment_matrix.csv", index_col=0)
meta = pd.read_csv(TABLE_DIR / "raw_peak_strain_treatment_metadata.csv")

meta.index = meta["sample"]
X_df["_res_sen"] = meta.loc[X_df.index, "label_res_sen"].values
X_df["_mech"]    = meta.loc[X_df.index, "mechanism"].values
X_df["_strain"]  = meta.loc[X_df.index, "strain_id"].values
X_df = X_df.dropna(subset=["_res_sen"])
print(f"  {X_df.shape[0]} rows matched")

feat_cols_all = [c for c in X_df.columns if not c.startswith("_")]
kruskal = pd.read_csv(TABLE_DIR / "raw_classic_biomarkers_kruskal_global.csv")
canon = set(kruskal["feature_chr"].astype(str))
feat_cols = [c for c in feat_cols_all if str(c) in canon]
if len(feat_cols) != len(canon):
    raise SystemExit(f"Features de R no presentes en la matriz: "
                     f"{len(canon)} esperadas, {len(feat_cols)} encontradas")
print(f"  Features (prevalencia maxima por clase >= 5%, via R): {len(feat_cols)}")
strain_ids = X_df["_strain"].to_numpy()

results = []

y_res = (X_df["_res_sen"] == "RES").astype(int).to_numpy()
df_res = run_conformal("res_sen", X_df[feat_cols].to_numpy().astype(float), y_res, strain_ids)
results.append(df_res)

mask = X_df["_mech"].isin(["KPC", "NDM"])
X_kpc = X_df[feat_cols][mask].to_numpy().astype(float)
y_kpc = (X_df.loc[mask, "_mech"] == "KPC").astype(int).to_numpy()
s_kpc = strain_ids[mask.to_numpy()]
df_kpc = run_conformal("kpc_ndm", X_kpc, y_kpc, s_kpc)
results.append(df_kpc)

out_df = pd.concat(results, ignore_index=True)
out_path = OUT_DIR / "conformal_prediction_raw_summary.csv"
out_df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print(out_df[["task", "balanced_accuracy", "auc", "coverage_95", "singleton_95"]].to_string(index=False))
