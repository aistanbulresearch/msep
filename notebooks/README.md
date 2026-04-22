# `msep` Reproducibility Notebooks

This folder groups notebooks by purpose. Every notebook is self-contained
(installs `msep`, fetches any required data, runs the analysis, saves the
figure). They are designed to run on **Google Colab Pro** (A100) without
local setup, but most will also work on a laptop with 8 GB RAM.

| Folder | Purpose |
|---|---|
| [`quickstart/`](quickstart/) | 5-minute smoke test for new users — uses the bundled synthetic dataset, no downloads required. |
| [`figures/`](figures/) | One notebook per paper figure (Figures 2–6). Each reproduces the four-panel layout from the bundled demo data and points to a real-data swap-in for paper-exact numbers. |
| [`validation/`](validation/) | Simulation-based ground-truth validation (Splatter / scDesign3) — WP-1. |
| [`pan_cancer/`](pan_cancer/) | CellxGene Census pan-cancer expansion — WP-3. |
| [`null_distributions/`](null_distributions/) | Gene-set permutation and cross-pathway null distributions — WP-5. |

## Opening a notebook in Colab

Click any of the "Open in Colab" badges in a notebook, or manually paste the
GitHub URL into Colab:

```
https://colab.research.google.com/github/aistanbulresearch/msep/blob/main/notebooks/<folder>/<notebook>.ipynb
```

## Full end-to-end example

For the complete chordoma analysis (all 106,584 cells, all figures in a
single walkthrough), see [`../examples/chordoma_msep.ipynb`](../examples/chordoma_msep.ipynb).
The notebooks in this directory are the **per-figure** and **per-work-package**
versions that the Genome Biology reviewer will exercise.
