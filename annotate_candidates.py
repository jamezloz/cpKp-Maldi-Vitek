#!/usr/bin/env python3
"""In silico UniProt annotation of the candidate masses (Table 3, Table S7) and decoy control.

One run of the procedure of Section 2.8: complete (non-fragment) K. pneumoniae entries,
[M+H]+ and [M+H]+ without the initial methionine, best match within +/-10 Da with
Swiss-Prot first. Candidate [M+2H]2+ peaks (Table S13) are also queried at the mass of
the singly charged species. Uniform and shifted decoys are processed identically.
Needs network access to the UniProt REST API. Writes to resultados/ANNOTATION/.
"""
from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from annotate_biomarkers_uniprot import best_match  # noqa: E402

BASE = Path(".")
OUT = BASE / "resultados/ANNOTATION"
CAND = BASE / "annotation/candidate_masses_tableS7.csv"
S13 = BASE / "resultados/PAPER_SUPP_TABLES/S13_double_charge_pairs.csv"
SEED = 42
GENERIC = ("KPHS", "KPN_", "NCTC", "KP1_", "SAMEA", "EL79", "BN49")


def row(mz: float, hit: dict | None) -> dict:
    if not hit:
        return {"mz": mz, "match": 0}
    return {"mz": mz, "match": 1, "hipotesis": hit["hipotesis"], "gen": hit["gen"],
            "proteina": hit["proteina"][:90], "accession": hit["accession"],
            "revisada": hit["revisada"], "masa_uniprot": hit["masa_uniprot"],
            "delta_da": hit["delta_da"],
            "gen_real": int(bool(hit["gen"]) and not hit["gen"].startswith(GENERIC))}


def resumen(rows: list[dict], conjunto: str) -> dict:
    n = len(rows)
    m = [r for r in rows if r["match"]]
    return {"conjunto": conjunto, "n": n, "con_match": len(m),
            "pct_match": round(100 * len(m) / n, 1),
            "pct_revisada": round(100 * sum(r["revisada"] == "si" for r in m) / n, 1),
            "pct_gen_nombrado": round(100 * sum(r["gen_real"] for r in m) / n, 1),
            "pct_delta_le_2": round(100 * sum(abs(float(r["delta_da"])) <= 2 for r in m) / n, 1)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    random.seed(SEED)
    cand = list(csv.DictReader(CAND.open()))
    masses = sorted({int(r["mz"]) for r in cand})
    pares = {round(float(r["mz_z2"])): float(r["mz_z1"]) for r in csv.DictReader(S13.open())}
    print(f"{len(cand)} entries, {len(masses)} distinct masses")

    reales = []
    for mz in masses:
        r = row(mz, best_match(mz))
        z1 = next((v for k, v in pares.items() if abs(k - mz) <= 2), None)
        if z1:
            h1 = best_match(z1)
            r.update({"z1_candidata": round(z1, 1),
                      "z1_gen": h1["gen"] if h1 else "", "z1_proteina": h1["proteina"][:90] if h1 else "",
                      "z1_accession": h1["accession"] if h1 else "", "z1_delta_da": h1["delta_da"] if h1 else "",
                      "z1_revisada": h1["revisada"] if h1 else ""})
        reales.append(r)
        print(f"  {mz}: {r.get('gen', '')} {r.get('proteina', 'no match')} {r.get('delta_da', '')}")

    lo, hi = min(masses), max(masses)
    unif = [row(m, best_match(m)) for m in (round(random.uniform(lo, hi)) for _ in masses)]
    desp = [row(m, best_match(m)) for m in
            (mz + random.choice((-1, 1)) * random.randint(150, 400) for mz in masses)]

    cols = ["mz", "match", "hipotesis", "gen", "proteina", "accession", "revisada", "masa_uniprot",
            "delta_da", "gen_real", "z1_candidata", "z1_gen", "z1_proteina", "z1_accession",
            "z1_delta_da", "z1_revisada"]
    for name, rows in (("annotation_candidates_v2.csv", reales),
                       ("annotation_decoys_uniform_v2.csv", unif),
                       ("annotation_decoys_shifted_v2.csv", desp)):
        with (OUT / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    res = [resumen(reales, "candidatas"), resumen(unif, "senuelo_uniforme"), resumen(desp, "senuelo_desplazado")]
    with (OUT / "annotation_decoy_summary_v2.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0]))
        w.writeheader()
        w.writerows(res)
    for r in res:
        print(r)


if __name__ == "__main__":
    main()
