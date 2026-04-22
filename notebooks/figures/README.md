# Per-Figure Reproducibility Notebooks

Each notebook in this folder reproduces one figure from Çavuş & Kuşkucu (2026).
All notebooks default to the bundled synthetic demo (`msep.datasets.load_example()`)
so a reviewer can go from `pip install msep` to a full figure in under 5 minutes
on a plain Colab CPU runtime. Each notebook also documents the real-data swap-in
for paper-exact numbers.

| Notebook | Reproduces | Quick-Colab |
|---|---|---|
| [`figure2_paradox.ipynb`](figure2_paradox.ipynb) | Multi-scale entropy paradox (entropy violin, CSC pathway CV, CV heatmap, paradox scatter) | [Open](https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/figures/figure2_paradox.ipynb) |
| [`figure3_pan_cancer.ipynb`](figure3_pan_cancer.ipynb) | Pan-cancer comparison of pathway CV (heatmap, EMT ranking, chordoma vs median) | [Open](https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/figures/figure3_pan_cancer.ipynb) |
| [`figure4_perturbation.ipynb`](figure4_perturbation.ipynb) | Cross-pathway pseudo-perturbation heatmap + gene-level CV for 15 defence genes | [Open](https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/figures/figure4_perturbation.ipynb) |
| [`figure5_xbp1.ipynb`](figure5_xbp1.ipynb) | XBP1 stress consolidation (distribution, CV by XBP1 group, pan-cancer scorecard) | [Open](https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/figures/figure5_xbp1.ipynb) |
| [`figure6_bulk_validation.ipynb`](figure6_bulk_validation.ipynb) | Cross-platform concordance (scRNA vs bulk, chordoma vs notochord, per-patient stability, gene-level Spearman) | [Open](https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/figures/figure6_bulk_validation.ipynb) |

## Running one

Colab CPU runtime is enough for every notebook here. Click the badge, then
**Runtime → Run all**. Each notebook writes the figure to the Colab working
directory as both PDF (vector, for the paper) and PNG (raster, for viewing).

## Demo vs. real data

The demo reproduces the paper's *qualitative* rankings on a 500-cell synthetic
AnnData: CSC highest per-cell entropy, EMT CV < ferroptosis CV < immune
evasion CV, XBP1-high consolidates all three defence pathways. Absolute
magnitudes are smaller than the 106,584-cell Arrieta 2025 cohort — the
demo dataset is designed for fast inspection of the MSEP output shape and
the pipeline wiring, not for biological publication.

For paper-exact numbers, each notebook's last markdown cell shows the
real-data substitution — typically one line that swaps
`msep.datasets.load_example()` for a loader that reads the Zenodo /
GEO / CellxGene Census cohort. The downstream code is unchanged.
