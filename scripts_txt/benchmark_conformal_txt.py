from __future__ import annotations

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

TABLE_DIR = Path("resultados/TXT_modular/tables")
OUT_DIR   = TABLE_DIR

def run_conformal(task: str, Xmat: np.ndarray, y: np.ndarray,
                  strain_ids: np.ndarray, n_outer: int = N_OUTER) -> pd.DataFrame:
    print(f"\n=== {task} ===")
    gss  = GroupShuffleSplit(n_splits=n_outer, test_size=0.2, random_state=RANDOM_STATE)
    rows = []

    for fold_i, (train_idx, test_idx) in enumerate(gss.split(Xmat, y, groups=strain_ids)):
        X_tr, X_te = Xmat[train_idx], Xmat[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
            continue

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

        y_pred = clf.predict(X_te)
        ba  = balanced_accuracy_score(y_te, y_pred)
        try:
            prob  = clf.predict_proba(X_te)
            pos_i = list(clf.classes_).index(1)
            auc   = roc_auc_score((y_te == 1).astype(int), prob[:, pos_i])
        except Exception:
            auc = np.nan

        fold_row = {"fold": fold_i, "balanced_accuracy": ba, "auc": auc}

        le       = LabelEncoder().fit(np.unique(y))
        y_te_enc = le.transform(y_te)
        results_conf = {}
        for alpha in ALPHAS:
            conf_level = 1 - alpha
            mapie = SplitConformalClassifier(
                clf, confidence_level=conf_level, prefit=True,
                random_state=RANDOM_STATE)
            mapie.conformalize(X_cal, y_cal)
            results_conf[conf_level] = mapie

        for alpha in ALPHAS:
            conf_level = 1 - alpha
            _, sets    = results_conf[conf_level].predict_set(X_te)
            sets_sq    = sets[:, :, 0]
            set_sizes  = sets_sq.sum(axis=1)
            covered    = np.array([bool(sets_sq[i, y_te_enc[i]]) for i in range(len(y_te))])
            lvl        = int(conf_level * 100)
            fold_row[f"coverage_{lvl}"]      = float(covered.mean())
            fold_row[f"singleton_{lvl}"]     = float((set_sizes == 1).mean())
            fold_row[f"empty_{lvl}"]         = float((set_sizes == 0).mean())
            fold_row[f"mean_set_size_{lvl}"] = float(set_sizes.mean())
        rows.append(fold_row)

    df = pd.DataFrame(rows)
    print(f"  BA={df['balanced_accuracy'].mean():.4f}  AUC={df['auc'].mean():.4f}")
    for alpha in ALPHAS:
        lvl = int((1 - alpha) * 100)
        print(f"  α={alpha}: coverage={df[f'coverage_{lvl}'].mean():.4f}  "
              f"singleton={df[f'singleton_{lvl}'].mean():.4f}")
    return (df.mean(numeric_only=True).to_frame().T
              .assign(task=task, model="rf", aggregation="strain_mean", n_folds=n_outer))

print("Loading TXT matrix and metadata...")
X_df = pd.read_csv(TABLE_DIR / "txt_matrix_aggregated.csv", index_col=0)
meta = pd.read_csv(TABLE_DIR / "txt_metadata_aggregated.csv")

print(f"  X shape: {X_df.shape}")
print(f"  meta shape: {meta.shape}")

kruskal = pd.read_csv(TABLE_DIR / "biomarkers_global_kruskal.csv")
canon = set(kruskal["feature"].astype(str))
feat_cols_f = [c for c in X_df.columns if str(c) in canon]
if len(feat_cols_f) != len(canon):
    raise SystemExit(f"Features de R no presentes en la matriz: "
                     f"{len(canon)} esperadas, {len(feat_cols_f)} encontradas")
print(f"  Features (prevalencia maxima por clase >= 5%, via R): {len(feat_cols_f)}")

Xmat = X_df[feat_cols_f].values.astype(float)

meta["strain_id"] = meta["sample"].str.split("__").str[0]
strain_ids = meta["strain_id"].values

results = []

mask_rs = meta["label_res_sen"].isin(["RES", "SEN"])
y_res   = (meta.loc[mask_rs, "label_res_sen"] == "RES").astype(int).values
s_res   = strain_ids[mask_rs.values]
df_res  = run_conformal("res_sen", Xmat[mask_rs.values], y_res, s_res)
results.append(df_res)

mask_kn = meta["mechanism"].isin(["KPC", "NDM"])
y_kn    = (meta.loc[mask_kn, "mechanism"] == "KPC").astype(int).values
s_kn    = strain_ids[mask_kn.values]
df_kn   = run_conformal("kpc_ndm", Xmat[mask_kn.values], y_kn, s_kn)
results.append(df_kn)

out_df   = pd.concat(results, ignore_index=True)
out_path = OUT_DIR / "conformal_prediction_txt_summary.csv"
out_df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print(out_df[["task", "balanced_accuracy", "auc", "coverage_95", "singleton_95"]].to_string(index=False))
