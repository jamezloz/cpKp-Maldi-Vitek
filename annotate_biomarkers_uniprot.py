from __future__ import annotations

import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = Path(".")
PANEL = BASE / "resultados/TXT_modular/tables/biomarker_panel_consolidado.csv"

PUBLICADAS = BASE / "resultados/referencias/MASTER_biomarker_annotations.csv"
OUT = BASE / "resultados/TXT_modular/tables/biomarker_panel_anotacion.csv"

TAXA = ("573", "72407")
WINDOW = 20
REPORT_TOL = 10
MET_LOSS = 131.04

RESIST_KW = ("beta-lactamase", "carbapenem", "efflux", "porin", "resistance",
             "multidrug", "antibiotic", "penicillin-binding", "quaternary ammonium")

def query(mass_lo: float, mass_hi: float) -> list[dict]:
    q = f"(taxonomy_id:{TAXA[0]} OR taxonomy_id:{TAXA[1]}) AND mass:[{int(mass_lo)} TO {int(mass_hi)}]"
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(q)
           + "&fields=accession,reviewed,gene_names,protein_name,mass,cc_subcellular_location"
           + "&format=tsv&size=50")
    try:
        with urllib.request.urlopen(url, timeout=40) as r:
            rows = r.read().decode("utf-8").splitlines()
    except Exception as e:
        print(f"    ! error de consulta: {e}", file=sys.stderr)
        return []
    if len(rows) < 2:
        return []
    return list(csv.DictReader(rows, delimiter="\t"))

def best_match(mz: float) -> dict | None:
    best = None
    for hypothesis, target in (("[M+H]+", mz - 1.008),
                               ("[M+H]+ sin Met", mz - 1.008 + MET_LOSS)):
        for rec in query(target - WINDOW, target + WINDOW):
            try:
                m = float(rec["Mass"].replace(",", ""))
            except (KeyError, ValueError):
                continue
            delta = m - target
            if abs(delta) > REPORT_TOL:
                continue
            reviewed = rec.get("Reviewed", "") == "reviewed"

            score = (0 if reviewed else 1, abs(delta))
            if best is None or score < best["_score"]:
                name = rec.get("Protein names", "")
                best = {
                    "_score": score,
                    "hipotesis": hypothesis,
                    "accession": rec.get("Entry", ""),
                    "revisada": "si" if reviewed else "no",
                    "gen": (rec.get("Gene Names", "") or "").split(" ")[0],
                    "proteina": name,
                    "masa_uniprot": round(m, 1),
                    "delta_da": round(delta, 1),
                    "localizacion": (rec.get("Subcellular location [CC]", "") or "")[:60],
                    "rel_resistencia": "si" if any(k in name.lower() for k in RESIST_KW) else "no",
                }
        time.sleep(0.3)
    return best

def main() -> None:
    panel = [r for r in csv.DictReader(PANEL.open()) if r["panel_principal"] == "TRUE"]
    ya = set()
    if PUBLICADAS.exists():
        for r in csv.DictReader(PUBLICADAS.open()):
            try:
                ya.add(float(r["mz"]))
            except (KeyError, ValueError):
                pass

    out = []
    for r in panel:
        mz = float(r["mz"])
        previa = any(abs(p - mz) <= 10 for p in ya)
        print(f"  {mz:8.0f} Da ({r['tarea']}, {r['nivel']})"
              + ("  [ya anotada en el paper]" if previa else ""))
        hit = best_match(mz)
        row = {"tarea": r["tarea"], "mz": int(mz), "nivel": r["nivel"],
               "anotada_previamente": "si" if previa else "no"}
        if hit:
            hit.pop("_score")
            row.update(hit)
            print(f"           -> {hit['gen'] or '(sin gen)'} | {hit['proteina'][:55]}"
                  f" | delta {hit['delta_da']} Da | {hit['hipotesis']}")
        else:
            row["proteina"] = "sin match a +/-10 Da"
            print("           -> sin match")
        out.append(row)

    cols = ["tarea", "mz", "nivel", "anotada_previamente", "gen", "proteina",
            "accession", "revisada", "masa_uniprot", "delta_da", "hipotesis",
            "localizacion", "rel_resistencia"]
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    con = sum(1 for r in out if r.get("accession"))
    print(f"\nanotadas {con} de {len(out)} masas del panel")
    print(f"guardado en: {OUT}")

if __name__ == "__main__":
    main()
