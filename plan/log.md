# MSEP — GB Submission Work Log

Each iteration appends a new section. Most recent iteration at the top. See
[`GB_SUBMISSION_WORKPLAN.md`](GB_SUBMISSION_WORKPLAN.md) for the overall plan
and decision log.

---

## Iteration 3 — 2026-04-22 · WP-2.3 per-figure reproducibility notebooks

**Duration:** ~1 session

### What landed

| Notebook | Purpose | Runs on demo? |
|---|---|---|
| `notebooks/figures/figure2_paradox.ipynb` | **Flagship** — 4-panel Figure 2 (entropy violin, CSC pathway CV bar, cell type × pathway heatmap, paradox scatter) with combined 2×2 publication layout | ✅ full |
| `notebooks/figures/figure3_pan_cancer.ipynb` | Figure 3 (pan-cancer heatmap, EMT ranking, chordoma vs median) — uses paper Table S5 values until WP-3 Census pipeline provides live matrix | ✅ full from paper values |
| `notebooks/figures/figure4_perturbation.ipynb` | Figure 4 panels A (pseudo-perturbation heatmap) + C (gene-level CV bars); panel B documented as illustrative schematic | ✅ full |
| `notebooks/figures/figure5_xbp1.ipynb` | Figure 5 (XBP1 distribution, pathway CV by XBP1 group with Δ annotations, pan-cancer consolidation scorecard) | ✅ full |
| `notebooks/figures/figure6_bulk_validation.ipynb` | Figure 6 (scRNA vs bulk ranking concordance, chordoma vs notochord, per-patient EMT-first, gene-level Spearman) | ✅ full |
| `notebooks/figures/README.md` | Index with Open-in-Colab links + demo-vs-real-data explanation | — |

### Pattern established

All notebooks follow the same structure: install cell → load demo data → `msep.profile(...)` → per-panel rendering → save PDF + PNG → real-data swap-in markdown block at the end. The `msep.profile()` call is the only stateful step; swapping in the Arrieta 2025 cohort changes nothing downstream because the plotting code is dataset-agnostic.

### Verification

- All 5 notebooks execute end-to-end locally (CPU, <5 min each) via `nbclient` smoke test — 0 errors
- Each panel writes both `.pdf` (vector) and `.png` (300 dpi raster) outputs
- Figure 4 perturbation schema aligned with actual `result.perturbation` columns (`delta_cv`, `perm_p`) after an initial column-name mismatch was caught by the smoke test

### Out of scope

- Panel B of Figure 4 (mutual-reinforcement schematic) is illustrative, not data-driven; documented as a vector-graphics task rather than generated from code.
- Real-data reproduction of Figures 3 and 5 panel C requires the WP-3 CellxGene Census pipeline; hard-coded paper values stand in until that lands.

---

## Iteration 2 — 2026-04-22 · PyPI release 1.1.0 + quickstart simplification

See [PR #4](https://github.com/aistanbulresearch/msep/pull/4) for details (merged to `main`).

---

## Iteration 1 — 2026-04-21 · WP-2.1 + WP-2.2 (bundled demo + Colab quickstart)

**Duration:** ~1 session · **Branch:** feature branch merged into `main` as PR #1

### What landed

| Area | Change |
|---|---|
| New module | `msep/datasets.py` — `load_example()` builds a deterministic 500-cell synthetic AnnData (4 cell types × 4 cancer-defense pathways + XBP1 + 120 background genes) that reproduces the paper's qualitative paradox |
| New tests | `tests/test_datasets.py` — 21 tests: structure, determinism, phenotype (CSC highest entropy, EMT CV < ferroptosis < immune in CSC, TBXT CV < 0.3 in CSC), XBP1 stress groups, end-to-end profile run |
| New notebook | `notebooks/quickstart/quickstart_colab.ipynb` — 8 code cells, executes end-to-end locally, Open-in-Colab badge, ≤5-min reviewer demo |
| New folder | `notebooks/{quickstart,figures,validation,pan_cancer,null_distributions}/` + `notebooks/README.md` — WP-aligned scaffold for future notebooks |
| Plan docs | `plan/GB_SUBMISSION_WORKPLAN.md` (consolidated work plan) + `plan/log.md` (this file) |
| README | Colab badge added; new "bundled example" quickstart block; new "Example datasets" API ref section; test badge 27/27 → 48/48 |
| Core fix | `MSEPResult.coordination_scores` property added — README advertised it but the attribute did not exist. Now wraps `scoring.coordination_score_table` and maps classification via `classify_paradox` |
| Exports | `msep.datasets` added to `msep/__init__.py` and `__all__` |
| Version | `1.0.1` → `1.1.0` in both `pyproject.toml` and `msep/__init__.py` |

### Decisions recorded

Stored in [`plan/GB_SUBMISSION_WORKPLAN.md §5`](GB_SUBMISSION_WORKPLAN.md#5-decisions-log-author-özge-2026-04-21) and in user memory (`memory/project_gb_submission.md`):

- ✅ **Pan-cancer scope locked: 8 extra cancer types** (HCC, PDAC, CRC, HNSCC, ovarian serous, mesothelioma, soft-tissue sarcoma, gastric adenocarcinoma) → total 20-type panel
- ✅ **Docs domain: GitHub Pages default** (`aistanbulresearch.github.io/msep`); no custom domain purchase
- ⏳ Deferred: BASiCS Option A vs B (remind at W4), Method Article reframing (remind at W5), Dual Shield reference (remind at WP-7)

### Verification

- `python -m pytest tests/` → **48 passed, 0 failed** (27 existing + 21 new)
- Notebook end-to-end local execution (minus the `pip install` cell) → **OK, no errors in any output cell**
- **Colab runtime verification (user-run, CPU, 2026-04-21)** — all cells executed to completion after install was pinned to the feature branch. XBP1 consolidation output reproduced the paper's pattern qualitatively: all three defense pathways show lower CV in XBP1-high vs XBP1-zero CSC (ferroptosis Δ=−0.064, immune_evasion Δ=−0.160, emt Δ=−0.049). Ranking preserved: `emt < ferroptosis < immune_evasion` within XBP1-zero CSC. Absolute magnitudes smaller than the full-cohort paper values as expected for 500-cell synthetic data.
- CI green on Python 3.9 / 3.10 / 3.11 / 3.12.
- Lint / style: module conforms to existing repo conventions (immutable specs via module-level dict, NB sampling, sparse CSR output, typed signatures).

### Files touched

**New:**
- `msep/datasets.py` (266 LOC)
- `tests/test_datasets.py` (240 LOC)
- `notebooks/quickstart/quickstart_colab.ipynb` (19 cells)
- `notebooks/README.md`
- `plan/GB_SUBMISSION_WORKPLAN.md`
- `plan/log.md` (this file)

**Modified:**
- `README.md` (+22 lines: Colab badge, quickstart block, datasets API, updated test badge)
- `msep/__init__.py` (version bump + datasets export)
- `msep/result.py` (+27 lines: `coordination_scores` property)
- `pyproject.toml` (version bump)

### Not yet done (reviewer-blocking items still pending)

1. Publish `msep==1.1.0` to PyPI (needs 2FA, user-side action)
2. Per-figure reproducibility notebooks (WP-2.3) — folder scaffold only
3. mkdocs docs site deployment (WP-2.5)
4. `nbmake` CI smoke test addition (WP-2.6)
5. Bundled `examples/chordoma_msep.clean.ipynb` (cleaned-output version of the existing 63 MB notebook)

### Next iteration candidates (ordered by GB impact)

- **WP-2.3** first figure notebook — `notebooks/figures/figure2_paradox.ipynb` using the bundled demo data (CPU-only, no Colab required)
- **WP-1 kickoff** — Splatter simulation scenarios A/B/C in `notebooks/validation/sim_splatter.ipynb` (needs Colab A100 later; can scaffold locally)
- **WP-3 kickoff** — CellxGene Census pull script in `notebooks/pan_cancer/fetch_and_profile.ipynb` for the 8 locked cancer types (needs Colab runtime)

---
