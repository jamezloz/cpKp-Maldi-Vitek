# cpKp-Maldi-Vitek

Analysis code for *From routine identification to resistance prediction: reproducible MALDI-TOF MS biomarkers and machine learning classify carbapenemase-producing* Klebsiella pneumoniae *across spectral data types* (Martínez-Lozano et al.).

Spectra (peak lists, mzML files, manifest and processed matrices) are deposited in Zenodo: https://doi.org/10.5281/zenodo.22738864.

## Layout

```
00_paper1_metadata.R        data path and file-name parsing (edit DATA_BASE here)
annotate_biomarkers_uniprot.py   UniProt in silico annotation (Table 3, Table S7)
annotate_decoy_control.py        decoy control for the mass-based annotation
scripts_txt/                peak-list (.txt) pipeline
scripts_raw/                mzML pipeline (MALDIquant, TopHat)
scripts_figures/            main-text and supplementary figures
scripts_supp/               CSV exports behind the supplementary tables
renv.lock, requirements.txt pinned R and Python environments
```

Scripts write to `resultados/TXT_modular/`, `resultados/RAW_modular/`, `resultados/PAPER_FIGURES/` and `resultados/PAPER_SUPP_TABLES/`, created on first run. Run everything from the repository root.

## Setup

1. Download and unzip the spectra from the Zenodo record, keeping the four folders `peaklists_siteA`, `peaklists_siteM`, `mzml_siteA`, `mzml_siteM`.
2. In `00_paper1_metadata.R` set `DATA_BASE` to the folder that contains those four folders. Strain, protocol, carbapenemase content, sequence type and replicate are parsed from the file names.
3. R: `renv::restore()` (R 4.5.2). Python: `pip install -r requirements.txt` (Python 3.13).

## Execution order

Peak-list dataset (`scripts_txt/`): `01` → `02` → `03`, `03b`, `03d` → `04` → `08` → `06c`, `06d`, `06e`, `06f` → `07d`, `07f` → `12` → `15` → `16` → `18` → `19` → `21` → `25` → `benchmark_conformal_txt.py`, `benchmark_dca_txt.py`.

RAW dataset (`scripts_raw/`): `17` → `01` → `02`, `02b` → `03` → `11` → `07` → `06d` → `24` → `23` → `13` → `12` → `15` → `26` → `28` → `benchmark_conformal_raw.py`, `benchmark_dca_raw.py`.

Then `annotate_biomarkers_uniprot.py` and `annotate_decoy_control.py` (both need network access to the UniProt REST API), the `scripts_supp/build_*.R` exports, and the figures: `scripts_figures/32_fig_pcoa_main.R` (Figure 1), `36_fig2_final_hybrid.R` (Figure 2), `31_fig_shap_composite.R` (Figure 3), `genera_sfig_suplementario.R [N]` (Figures S1–S8).

Random seeds are fixed inside the scripts; the bootstrap and cross-validation loops (100–200 iterations) take from minutes to a few hours per script.

## Licence

MIT. Please cite the paper when using this code.
