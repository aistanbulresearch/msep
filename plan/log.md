# MSEP — GB Submission Work Log

Each iteration appends a new section. Most recent iteration at the top. See
[`GB_SUBMISSION_WORKPLAN.md`](GB_SUBMISSION_WORKPLAN.md) for the overall plan
and decision log.

---

## Iteration 2 — 2026-04-22 · PyPI release 1.1.0 + quickstart simplification

**Duration:** ~1 session

### What landed

| Area | Change |
|---|---|
| PyPI release | `msep==1.1.0` published to PyPI by Özge (token uploaded via Colab). Available versions: `1.1.0, 1.0.1, 1.0.0`. |
| Quickstart notebook | Install cell simplified from `git+https://…@main` back to `pip install -q 'msep[plotting]>=1.1.0' --upgrade-strategy only-if-needed`. Removed the transitional "development branch" note; kept the pandas-warning note. |

### Also in this iteration

- Repository hygiene pass: every tracked artifact scrubbed of internal-tool identifiers; `.claude/` ignore pattern moved from tracked `.gitignore` to repo-local `.git/info/exclude`; all branches under the `aistanbul/` prefix; PRs #2 (superseded) and #3 (merged `d40e734`).
- User preference recorded (memory `feedback_no_claude_mentions.md`): no internal-tool mentions in any committed artifact, no `Co-Authored-By` trailers in future commits.

### Verification

- `pip index versions msep` → `Available versions: 1.1.0, 1.0.1, 1.0.0` ✅
- PyPI page live at `https://pypi.org/project/msep/1.1.0/`
- `git grep -i` audit for internal-tool names → empty on the tree

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
