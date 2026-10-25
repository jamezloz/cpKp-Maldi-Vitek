import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "PATH_TO_SPECTRA")
OUT = Path("resultados/CONFOUNDER")
OUT.mkdir(parents=True, exist_ok=True)

RES = {"KPC", "NDM", "KPC-NDM", "KPC-OXA48", "NDM-OXA48"}
HDR = re.compile(r"# Data:\s*(IVD|RUO)_\w+?_pl-\w+?_(\d+\.\d+\.\d+)_\S+\s+(\d{4}-\d{2}-\d{2})")

man = pd.read_csv(DATA / "spectra_manifest.csv", dtype={"strain": str})


def peaklist_table() -> pd.DataFrame:
    rows = []
    for f in man[man.data_type == "peak_list"].itertuples():
        m = HDR.search((DATA / f.folder / f.file).read_text(errors="ignore"))
        rows.append(dict(file=f.file, folder=f.folder, site=f.site, strain=f.strain,
                         mech=f.mechanism, mode=m.group(1), version=m.group(2), month=m.group(3)[:7]))
    d = pd.DataFrame(rows)
    d["cls"] = np.where(d.mech == "SEN", "CSKp", "CRKp")
    return d


def mzml_table() -> pd.DataFrame:
    rows = []
    for f in man[(man.data_type == "mzML") & (man.site == "A") & man.mechanism.isin(RES | {"SEN"})].itertuples():
        h = (DATA / f.folder / f.file).open(errors="ignore").read(9000)
        laser = re.search(r'value="([\d.]+)" name="laser energy"', h)
        date = re.search(r'value="(\d{6})\d*"[^>]*name="acquisitionGroup preparation date"', h)
        rows.append(dict(file=f.file, strain=f.strain, mech=f.mechanism,
                         laser=float(laser.group(1)) if laser else np.nan,
                         month=int(date.group(1)) if date else np.nan))
    d = pd.DataFrame(rows).dropna()
    d["cls"] = np.where(d.mech == "SEN", "CSKp", "CRKp")
    return d


def single_predictor_auc(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for task, col, pos, neg in [("CRKp/CSKp", "cls", "CRKp", "CSKp"), ("KPC/NDM", "mech", "KPC", "NDM")]:
        s = d[d[col].isin([pos, neg])]
        y = (s[col] == pos).astype(int)
        rows.append(dict(task=task, spectra=len(s),
                         auc_laser=round(max(a := roc_auc_score(y, s.laser), 1 - a), 3),
                         auc_month=round(max(b := roc_auc_score(y, s.month), 1 - b), 3)))
    return pd.DataFrame(rows)


def presence_matrix(d: pd.DataFrame) -> np.ndarray:
    edges = np.exp(np.arange(np.log(2000), np.log(20000), np.log(1.001)))
    X = np.zeros((len(d), len(edges) - 1), np.uint8)
    for i, f in enumerate(d.itertuples()):
        mz = np.loadtxt(DATA / f.folder / f.file, comments="#", usecols=0, ndmin=1)
        j = np.searchsorted(edges, mz) - 1
        X[i, j[(j >= 0) & (j < X.shape[1])]] = 1
    return X[:, X.mean(0) >= 0.02]


def blocked_auc(X, y, groups) -> float:
    k = min(5, len(np.unique(groups[y == 1])), len(np.unique(groups[y == 0])))
    pred = np.zeros(len(y))
    for tr, te in StratifiedGroupKFold(k, shuffle=True, random_state=0).split(X, y, groups):
        m = RandomForestClassifier(300, class_weight="balanced", random_state=0, n_jobs=-1).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, pred)


def single_factor_contrasts(d: pd.DataFrame) -> pd.DataFrame:
    X = presence_matrix(d)
    A, M = d.site == "A", d.site == "M"
    old = (d["mode"] == "IVD") & (d.version == "1.0.8")
    shared = d.strain.isin(set(d[A].strain) & set(d[M].strain))
    designs = [
        ("carbapenemase type, same months and configuration", "KPC vs NDM, Site A, Dec 2024-Feb 2025",
         A & old & d.mech.isin(["KPC", "NDM"]) & d.month.isin(["2024-12", "2025-01", "2025-02"]), lambda s: s.mech == "KPC"),
        ("acquisition month, same class and configuration", "KPC, Site A, Sep 2024 vs Feb 2025",
         A & old & (d.mech == "KPC") & d.month.isin(["2024-09", "2025-02"]), lambda s: s.month == "2025-02"),
        ("acquisition month, same class and configuration", "NDM, Site A, Jul-Aug 2024 vs Dec 2024-Feb 2025",
         A & old & (d.mech == "NDM") & d.month.isin(["2024-07", "2024-08", "2024-12", "2025-01", "2025-02"]),
         lambda s: s.month.isin(["2024-12", "2025-01", "2025-02"])),
        ("acquisition month, same class and configuration", "CSKp, Site M, Nov 2024-Apr 2025 vs Oct 2025",
         M & (d.cls == "CSKp"), lambda s: s.month == "2025-10"),
        ("site, same strains", "shared CRKp, Site A vs Site M (same configuration)",
         shared & (d.cls == "CRKp") & old, lambda s: s.site == "A"),
        ("site, same strains", "shared CSKp, Site A vs Site M",
         shared & (d.cls == "CSKp"), lambda s: s.site == "A"),
    ]
    rows = []
    for factor, contrast, mask, label in designs:
        mask = mask.to_numpy()
        s = d[mask]
        y = label(s).astype(int).to_numpy()
        rows.append(dict(factor=factor, contrast=contrast, spectra=int(mask.sum()),
                         strains=s.strain.nunique(), auc=round(blocked_auc(X[mask], y, s.strain.to_numpy()), 3)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pl = peaklist_table()
    design = pd.crosstab([pl.site, pl["mode"], pl.version, pl.month], pl.cls)
    design.to_csv(OUT / "acquisition_design.csv")
    print(design.to_string(), "\n")
    t1 = single_predictor_auc(mzml_table())
    t1.to_csv(OUT / "single_predictor_auc.csv", index=False)
    print(t1.to_string(index=False), "\n")
    t2 = single_factor_contrasts(pl)
    t2.to_csv(OUT / "single_factor_contrasts.csv", index=False)
    print(t2.to_string(index=False))
