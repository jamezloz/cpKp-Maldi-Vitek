from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from annotate_biomarkers_uniprot import best_match

BASE = Path(".")
ANOT = BASE / "resultados/TXT_modular/tables/biomarker_panel_anotacion.csv"
OUT_DET = BASE / "resultados/TXT_modular/tables/annotation_decoy_control_detail.csv"
OUT_SUM = BASE / "resultados/TXT_modular/tables/annotation_decoy_control_summary.csv"

SEED = 42
GENERIC_PREFIX = ("KPHS", "KPN_", "NCTC", "KP1_", "SAMEA", "EL79", "BN49")


def is_named(gene: str) -> bool:
    return bool(gene) and not gene.startswith(GENERIC_PREFIX)


def classify(hit: dict | None) -> dict:
    if not hit:
        return {"match": 0, "revisada": 0, "gen_real": 0, "delta2": 0}
    return {
        "match": 1,
        "revisada": int(hit["revisada"] == "si"),
        "gen_real": int(is_named(hit["gen"])),
        "delta2": int(abs(float(hit["delta_da"])) <= 2),
    }


def run_set(masses: list[float], label: str) -> list[dict]:
    rows = []
    for i, mz in enumerate(masses, 1):
        hit = best_match(float(mz))
        r = {"conjunto": label, "mz": int(mz)}
        r.update(classify(hit))
        if hit:
            r.update({"gen": hit["gen"], "proteina": hit["proteina"][:70],
                      "delta_da": hit["delta_da"], "revisada_txt": hit["revisada"]})
        rows.append(r)
        if i % 10 == 0:
            print(f"    {label}: {i}/{len(masses)}")
    return rows


def main() -> None:
    random.seed(SEED)
    real = sorted({int(float(r["mz"])) for r in csv.DictReader(ANOT.open())})
    lo, hi = min(real), max(real)
    print(f"panel real: {len(real)} masas, rango {lo}-{hi} Da\n")

    uniform = [round(random.uniform(lo, hi)) for _ in real]
    shifted = [m + random.choice([-1, 1]) * random.randint(150, 400) for m in real]
    shifted = [max(lo, min(hi, m)) for m in shifted]

    rows = []
    rows += run_set(real, "panel_real")
    rows += run_set(uniform, "senuelo_uniforme")
    rows += run_set(shifted, "senuelo_desplazado")

    with OUT_DET.open("w", newline="") as fh:
        cols = ["conjunto", "mz", "match", "revisada", "gen_real", "delta2",
                "gen", "proteina", "delta_da", "revisada_txt"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    summary = []
    for label in ("panel_real", "senuelo_uniforme", "senuelo_desplazado"):
        s = [r for r in rows if r["conjunto"] == label]
        n = len(s)
        summary.append({
            "conjunto": label, "n_masas": n,
            "con_match_pct": round(100 * sum(r["match"] for r in s) / n, 1),
            "revisada_pct": round(100 * sum(r["revisada"] for r in s) / n, 1),
            "gen_real_pct": round(100 * sum(r["gen_real"] for r in s) / n, 1),
            "delta_leq_2Da_pct": round(100 * sum(r["delta2"] for r in s) / n, 1),
        })
    with OUT_SUM.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)

    print("\n=== CONTROL SENUELO ===")
    hdr = f"{'conjunto':22s} {'n':>4s} {'match%':>8s} {'revisada%':>10s} {'gen real%':>10s} {'delta<=2%':>10s}"
    print(hdr); print("-" * len(hdr))
    for r in summary:
        print(f"{r['conjunto']:22s} {r['n_masas']:4d} {r['con_match_pct']:8.1f} "
              f"{r['revisada_pct']:10.1f} {r['gen_real_pct']:10.1f} {r['delta_leq_2Da_pct']:10.1f}")

    pr = summary[0]
    peor = max(summary[1:], key=lambda r: r["con_match_pct"])
    print()
    if pr["con_match_pct"] - peor["con_match_pct"] < 10:
        print("LECTURA: los senuelos casan al mismo ritmo que el panel. La coincidencia")
        print("de masa NO identifica proteina en este proteoma; las asignaciones son")
        print("candidatos compatibles y deben presentarse como tales.")
    else:
        print("LECTURA: el panel casa por encima del fondo; la anotacion aporta senal.")
    print(f"\nguardado en:\n  {OUT_SUM}\n  {OUT_DET}")


if __name__ == "__main__":
    main()
