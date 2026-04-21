# MSEP — Genome Biology Submission Work Plan

**Target journal:** Genome Biology (Method Article)
**Manuscript:** `Kordoma_Entropi_Makale (7).pdf` (26 pp, v7, 21 Apr 2026)
**Package:** `msep` v1.0.1 on PyPI (27/27 tests, CI green)
**Timeline:** 4–6 weeks to submission
**Owner:** Özge A. Çavuş (lead), Ayşegül Kuşkucu (co-author)

---

## 0. Current State Audit (as of 2026-04-21)

### What the paper already has
| Item | Status |
|---|---|
| Primary dataset (6 chordomas, 106,584 cells, Arrieta 2025) | ✅ complete |
| Bulk validation (GSE205457, Halvorsen 2023, 6 tumors + 2 notochord) | ✅ complete |
| Pan-cancer comparison (12 types, 10x + Smart-seq2) | ✅ complete |
| scVI Bayesian variance decomposition (BDR > 0.96 all pathways) | ✅ complete |
| Pseudo-perturbation (6 shield genes, 500 permutations) | ✅ complete |
| XBP1 consolidation (9 cancer types, 5/9 full consolidation) | ✅ complete |
| Fano-factor supplementary control | ✅ complete |
| Data & code availability statements | ✅ complete |

### What the repo already has
| Item | Status |
|---|---|
| `msep` v1.0.1 on PyPI, MIT license | ✅ |
| Module layout: core/entropy/coordination/perturbation/bayesian/scoring/pathways/plotting/result | ✅ (2,752 LOC, none over 440) |
| Unit tests: 27/27 passing | ✅ |
| GitHub Actions CI (Python 3.9–3.12) | ✅ |
| Example notebook: `examples/chordoma_msep.ipynb` | ✅ (single notebook; not per-figure) |
| MSigDB integration (12 collections, 33k+ sets) | ✅ |
| Coordination Score (single-number summary) | ✅ |
| ReadTheDocs / mkdocs site | ⚠️ **missing** (declared in URLs, not deployed) |
| Quickstart Colab notebook with badge | ⚠️ **missing** |
| Per-figure reproducibility notebooks | ⚠️ **missing** |
| Bundled example AnnData (`msep.datasets.load_example()`) | ⚠️ **missing** |
| Vignette / extended tutorial | ⚠️ **missing** |

### What the report (`gb_hazirlama_raporu.md`) flags as GB-blockers
1. **P0 — Simulation validation** (Splatter / scDesign3 ground truth) — *missing*
2. **P0 — Software hardening** (quickstart, vignette, per-figure repro, bundled demo data, docs site)
3. **P1 — Pan-cancer expansion** (+5–8 cancer types via CellxGene Census → 17–20 total)
4. **P1 — BASiCS comparison** (discussion text minimum; optional residual-overdispersion analysis)
5. **P2 — Gene-set-permutation null distribution** (1,000 perms, z-score vs. random gene sets)
6. **P2 — Method Article reframing** (title, abstract, results order)
7. **P3 — Dual Shield reference decision** (cut or reframe)

---

## 1. Priority Matrix

| # | Work package | GB impact | Effort | Duration | Priority | Owner |
|---|---|---|---|---|---|---|
| WP-1 | Splatter/scDesign3 simulation validation → Fig S12/S13 | Critical | Medium | 1–2 wk | **P0** | Özge |
| WP-2 | `msep` package hardening (Colab quickstart, docs site, per-figure notebooks, bundled data, vignette) | Critical | Medium | 1 wk | **P0** | Özge |
| WP-3 | Pan-cancer expansion (+5–8 types via CellxGene Census) → update Fig 3, 5C | High | Low | 3–5 d | **P1** | Özge |
| WP-4 | BASiCS comparison (Discussion write-up + optional residual overdispersion on CSC) | High | Medium | 1 wk | **P1** | Özge |
| WP-5 | Gene-set permutation null (1,000 random gene-matched sets) → Supplementary Fig | Medium | Low | 2–3 d | P2 | Özge |
| WP-6 | Method Article reframing (title, abstract, results order, introduction gap section) | Medium | Medium | 3–5 d | P2 | Özge + Ayşegül |
| WP-7 | Dual Shield citation decision (prefer removal) + TBXT-excluded EMT CV robustness | Low | Minimal | 1 d | P3 | Özge |
| WP-8 | Cover letter, BMC Minimum Standards Checklist, ORCID sweep, APC waiver drafting | Submission | Low | 2 d | P1 | Özge |

---

## 2. Work Package Details

### WP-1 — Simulation Validation (P0)

**Goal:** Prove MSEP detects ground-truth "individually diverse, collectively disciplined" phenotype under controlled synthetic conditions. Rebuts reviewer concern that low EMT CV is a technical artifact.

**Deliverables**
- `notebooks/validation/sim_splatter.ipynb` — 3 scenarios × 50k cells × 4 pathways × 3 seeds
- `src/validation/simulate.py` — parameterised Splatter runner
- Supplementary Fig **S12** (phenotype recovery under ground truth), **S13** (CV–mean decoupling)
- Methods §2.x subsection + Results one-paragraph + Discussion half-paragraph

**Scenarios (must match paper's reported CSC profile)**
| Scenario | Per-cell entropy | EMT CV | Immune CV | Label |
|---|---|---|---|---|
| A (target) | high (≥9.5 bits) | low (≤5) | high (≥8) | "individually diverse, collectively disciplined" |
| B (negative control) | low | low | low | "uniform and constrained" |
| C (negative control) | high | high | high | "chaotic" |

**Implementation notes**
- Use `splatter` via `rpy2` (CRAN `splatter::splatSimulate`) — more mature than `splatter-py`
- `delta` (dispersion) knob per pathway: EMT δ=0.3, Immune δ=3.0, Ferroptosis δ=1.0, HK δ=0.2
- `dropout.type='experiment'` with capture efficiency matched to Arrieta 10x 3′ v3 data
- 6 batches (patient proxy), `batch.facLoc=0.05` to mimic Harmony-corrected residual batch
- Primary metric: recover MSEP ranking (EMT < ferroptosis < immune CV) with > 95% accuracy across seeds
- Report `true_delta` vs. `recovered_pathway_CV` scatter → proof of calibration

**Compute:** Colab Pro A100, ~3–4 h total (3 scenarios × 3 seeds × ~20 min Splatter + ~10 min MSEP).

**Acceptance gate**
- MSEP recovers scenario labels blindly from CV profile alone → include as Supplementary Fig legend
- scVI-corrected CV vs. raw CV scatter on simulated data shows no spurious divergence (Fig S13)

---

### WP-2 — `msep` Package Hardening (P0)

**Goal:** A reviewer must go from `pip install msep` → reproducing Figure 2 in **≤30 minutes**.

**Sub-tasks**

| # | Task | Acceptance |
|---|---|---|
| 2.1 | Bundle example AnnData (`msep.datasets.load_example()`) | ≈500 cells, 4 pathways, raw_counts layer, patient_id, cell_type; ≤2 MB |
| 2.2 | Colab quickstart notebook with badge | Open-in-Colab button in README, installs msep, loads example, shows paradox + CV heatmap, total runtime ≤5 min |
| 2.3 | Per-figure reproducibility notebooks | `notebooks/figure2.ipynb` … `notebooks/figure6.ipynb` each self-contained (install + data-download + compute + save figure) |
| 2.4 | Vignette (advanced) | `docs/vignette.md` — custom pathway definition, BYO-AnnData walkthrough, XBP1 extension |
| 2.5 | mkdocs-material docs site | Deploy to `aistanbulresearch.github.io/msep` via GitHub Actions; every public function has docstring + example |
| 2.6 | CI extensions | Add `pip install msep` clean-install job, notebook execution smoke test via `nbmake` |
| 2.7 | CHANGELOG + version bump | v1.1.0 "GB submission release" |

**Risk mitigation**
- Test the full "cold-start reviewer" path on a fresh Colab Pro runtime before each submission-draft iteration
- Add "Reproduce the paper" section to README linking the 5 figure notebooks + Zenodo data DOI

---

### WP-3 — Pan-Cancer Expansion (P1)

**Goal:** Extend from 12 → 17–20 cancer types to pre-empt "is your framework actually generalisable?" reviewer push.

**Additional types (CellxGene Census v2025-11-08)** — **locked at 8 types (broadest scope)**, total panel 20 cancer types post-expansion
1. Hepatocellular carcinoma (HCC)
2. Pancreatic ductal adenocarcinoma (PDAC)
3. Colorectal carcinoma (CRC)
4. Head & neck squamous cell carcinoma (HNSCC)
5. Ovarian serous carcinoma
6. Mesothelioma
7. Soft tissue sarcoma
8. Gastric adenocarcinoma

**Pipeline (per cancer type)**
1. Census query, `is_primary_data==True`, 10x-only
2. Malignant filter (author-labeled or `cell_type ∈ malignant_lexicon`)
3. Subsample to ≤8,000 cells (match existing protocol, seed=42)
4. `msep.profile(adata, pathways="cancer_defense", cell_type_key=<malignant_label>)`
5. Append to Fig 3 heatmap + pan-cancer median tables + XBP1 consolidation scorecard

**Deliverables**
- `notebooks/pan_cancer_extension.ipynb`
- Updated Fig 3 (17–20 rows), Fig 5C scorecard
- Updated Supp Table S5 (pan-cancer CV complete table)
- Results §3.4 reword: "Among 17 cancer types spanning 8 CellxGene Census queries and 3 GEO datasets, chordoma CSC exhibited the lowest EMT CV…"

**Compute:** ~3–5 h Colab Pro (Census query bottleneck, not MSEP).

---

### WP-4 — BASiCS Comparison (P1)

**Decision tree** — lead with **Option A (discussion text)**, upgrade to **Option B (analysis)** only if a time budget remains after WP-1 & WP-3.

**Option A (discussion, 1 day)**
Add a Discussion paragraph contrasting MSEP + scVI with BASiCS:
- scVI models batch, library size, dropout jointly with a deep generative prior; BDR > 0.96 already addresses the technical/biological decomposition BASiCS was designed for.
- BASiCS works at the gene level; MSEP's contribution is the pathway-level aggregation and cross-pathway coordination, which BASiCS does not natively provide.
- Cite Lin 2023 (GB benchmark on raw-CV in 10x data) and Eling 2018 (BASiCS).

**Option B (analysis, 1 week, rpy2 + BASiCS R)**
- Run `BASiCS_MCMC` on the 6,730 CSC_TBXT+ cells (CPU-bound, 4–8 h)
- Compute residual over-dispersion per gene
- Recompute pathway-level ranking using residual OD → confirm EMT < ferroptosis < immune_evasion preserved
- New Supp Fig: scVI BDR ordering vs. BASiCS residual-OD ordering (Spearman)

**Acceptance gate for Option B:** if BASiCS ordering diverges from scVI + raw CV, **do not include** — Option A becomes the final submission content. (A divergent result opens a question we cannot resolve inside the 6-week window.)

---

### WP-5 — Gene-Set Permutation Null (P2)

**Goal:** Directly answer "what if pathway assignment were random?" — complement the existing 500-iter pseudo-perturbation with a 1,000-iter gene-set null.

**Two null distributions**
1. **Gene-set permutation:** replace each pathway's genes with N mean-expression-matched random genes; recompute CV. 1,000 permutations × 3 pathways. Report z-score for observed pathway CV vs. null.
2. **Cross-pathway dependency null:** shuffle gene-to-pathway assignments across the three defense sets; test whether the VIM-high → ferroptosis-CV-drop effect survives.

**Deliverables**
- `notebooks/null_distributions.ipynb`
- Supp Fig with null histograms + observed-value z-scores
- One Methods paragraph + one Results sentence in §3.5

**Compute:** <1 h on Colab (vectorised numpy on sparse matrix).

---

### WP-6 — Method Article Reframing (P2)

**Goal:** Shift the submission genre so chordoma biology becomes a *compelling case study* inside a methods-led narrative — matches GB's method-article expectations (BASiCS, sctransform, Scrublet all followed this template).

**Structural edits**

| Element | Current | Proposed |
|---|---|---|
| Title | "Multi-Scale Entropy Profiling Reveals Individually Diverse but Collectively Disciplined Transcriptional Programs in Chordoma Stem Cells" | "Multi-scale entropy profiling: a framework for quantifying pathway-selective population coordination in single-cell transcriptomics" |
| Abstract §1 | Chordoma biology lead | Framework lead ("Existing entropy analyses operate at a single scale…"); chordoma relegated to "We demonstrate this framework on 106,584 cells from six chordomas and generalise to 17 cancer types" |
| Introduction | Chordoma-first | Gap-first: SCENT / Kinker meta-programs / Patel constraint model → limitations → MSEP |
| Results order | (1) landscape, (2) paradox, (3) pan-cancer, (4) perturbation, (5) XBP1, (6) bulk | (1) Framework + simulation validation, (2) Chordoma case study, (3) Pan-cancer generalisation, (4) XBP1 consolidation |
| Discussion | Chordoma-heavy | Framework-vs-framework (BASiCS, scVI, SCENT, Kinker) first; chordoma-specific second |

**Keep** the "individually diverse, collectively disciplined" phrase — it is the signature finding and survives reframing.

---

### WP-7 — Minor Revisions (P3)

1. **Dual Shield reference (Section 4.3 line 570–571)**
   - Recommended action: **remove the sentence** (paper stands without it; ESHG poster is not peer-reviewed; keeps submission tidy)
2. **TBXT-exclusion EMT CV robustness**
   - Recompute EMT pathway CV with TBXT excluded → confirm ranking preserved in all 12 (17 after WP-3) cancer types
   - Add one sentence + Supp Figure to Limitations §4.5
3. **Notochord n=2 framing**
   - Rephrase §3.7 last paragraph as "hypothesis-generating"; optionally mine GTEx developmental datasets for additional notochord-adjacent tissue (low priority)

---

### WP-8 — Submission Packaging (P1)

**Checklist**

| Item | Status | Action |
|---|---|---|
| Data public (Zenodo + GEO) | ✅ | — |
| Code public (GitHub `msep`) | ✅ | — |
| OSI-compliant license (MIT) | ✅ | — |
| Software installs & runs in ≤30 min | 🟡 | Delivered by WP-2 |
| Per-figure reproducibility scripts | 🟡 | Delivered by WP-2.3 |
| Synthetic validation | 🟡 | Delivered by WP-1 |
| Cover letter | ❌ | Draft (WP-8) |
| Data Availability Statement | ✅ | — |
| ORCID IDs both authors | ❓ | Verify |
| BMC Minimum Standards Reporting Checklist | ❌ | Fill during WP-8 |
| APC waiver application | ❌ | Draft during WP-8 (Genome Biology discretionary waiver) |

**Cover letter skeleton** (~300 words, 3 paragraphs per report §10):
- Hook: single-scale entropy limitation → MSEP fills it
- Significance: 12+ (17 post-WP-3) cancer types, bulk validation, open-source tool
- Why GB: method-article fit (BASiCS, sctransform precedent), reproducibility tier met

---

## 3. Timeline (6-week plan)

| Week | Parallelisable tasks | Serial blockers | Outputs |
|---|---|---|---|
| **W1** | WP-1 Splatter setup · WP-2.1/2.2 bundled data + Colab notebook · WP-7.2 TBXT-excluded CV | — | Sim scenario A draft; Colab quickstart live; TBXT robustness figure |
| **W2** | WP-1 finish 3 scenarios · WP-2.3 per-figure notebooks · WP-3 CellxGene pulls | WP-1 Scenario A must land before S13 | Supp S12/S13 committed; notebooks 2–6 runnable; 5 new cancer types processed |
| **W3** | WP-3 integrate → Fig 3/5C · WP-2.5 mkdocs deploy · WP-5 null distribution | — | Updated Fig 3 (17 rows); docs site live; null histograms figure |
| **W4** | WP-4 BASiCS discussion draft (+ optional analysis) · WP-2.4 vignette · WP-2.6 CI smoke tests | WP-4 Option A/B decision gate end of W4 | Revised Discussion; `msep` v1.1.0 on PyPI |
| **W5** | WP-6 Method-Article reframing · WP-7.1 Dual Shield removal · WP-7.3 notochord phrasing | Manuscript rewrite depends on W1–W4 outputs | Revised manuscript (v8) |
| **W6** | WP-8 cover letter · BMC checklist · ORCID + APC waiver · final read-through | — | **Submission** |

**Parallelism:** WP-1, WP-2, WP-3, WP-5 can all run in parallel across W1–W3. WP-4, WP-6 must come after the primary-content packages.

---

## 4. Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Splatter simulation fails to reproduce MSEP ranking | Low | High | Scenario A parameter sweep in W1; fall back to scDesign3 if needed |
| CellxGene Census query rate-limited or 10x-only filter too strict | Medium | Low | Pre-cache queries locally; document fallback to 3 GEO datasets |
| BASiCS Option B ranking diverges from scVI | Medium | Medium | Acceptance gate — drop Option B, keep Option A only |
| Colab notebook flakiness on fresh runtime | Medium | Medium | Pin all dependency versions; add `nbmake` CI smoke test |
| Method-Article reframing delays writing to W5 and compresses submission | Medium | High | Start abstract + introduction skeleton in W3 alongside WP-3 data work |
| APC (~$5,490) waiver denied | Medium | High | Draft waiver argument early (W3), prepare fallback co-funding source |

---

## 5. Decisions Log (author: Özge, 2026-04-21)

| # | Question | Decision | Status |
|---|---|---|---|
| 1 | BASiCS Option A (discussion) vs. Option B (analysis) | **Defer** — evaluate when WP-4 week reached | ⏳ PENDING — remind before W4 start |
| 2 | Dual Shield reference (Section 4.3) — keep / reframe / remove | **Defer** — decide later | ⏳ PENDING — remind during WP-7 |
| 3 | Method Article reframing (title, abstract, results order) | **Defer** — decide later | ⏳ PENDING — remind before W5 start |
| 4 | Pan-cancer expansion scope | **Broadest — 8 additional cancer types** | ✅ LOCKED (WP-3) |
| 5 | Docs site domain | No custom domain for now; purchase only if needed | ✅ LOCKED — use GitHub Pages default for WP-2.5 |

### Reminder triggers

Claude will proactively surface the pending decisions at these moments:
- **Start of Week 4** (WP-4 kick-off) → remind about #1 (BASiCS Option A vs B)
- **Start of Week 5** (WP-6 kick-off) → remind about #3 (Method Article reframing)
- **During WP-7 Minor Revisions** → remind about #2 (Dual Shield reference)

### Implications for the plan

- **WP-3 is locked to 8 extra types:** HCC, PDAC, CRC, HNSCC, ovarian serous, mesothelioma, soft-tissue sarcoma, gastric adenocarcinoma → total pan-cancer panel **20 types** (12 existing + 8 new). Compute budget accordingly: ~5–6 h Colab Pro for Census pulls.
- **WP-2.5 docs deployment target:** `https://aistanbulresearch.github.io/msep/` (GitHub Pages default). No DNS / custom-domain work needed in the 6-week window.
- **WP-4, WP-6, WP-7 hold partial deferrals:** structure the week calendar so these packages begin only after decisions are confirmed; provide fallback default if author is unavailable at reminder time (default = report recommendation: Option A for #1, remove for #2, reframe for #3).

---

## 6. Tracking

- Progress tracker: update this file weekly; move completed items to a `## Done` section at the bottom
- GitHub Issues: open one issue per WP, labelled `gb-submission`, with checkboxes for sub-tasks
- Branch strategy: `feature/gb-wp1-sim-validation`, `feature/gb-wp2-docs`, etc.; rebase onto `main` after review
- Commit prefix: `paper(gb):` for manuscript-side edits, standard conventional commits for code
