# Manuscript revision notes

Reproducibility findings from rebuilding the figures and supplementary materials of
Çavuş & Kuşkucu (2026) against the analysis output in `processed/`. Each item is a
change the manuscript text or a figure needs so that it matches the data behind it.

Every numeric claim listed here is now guarded by a test in `tests/`, so once the
text is updated it will not silently drift again. Nothing in this document changes a
scientific conclusion — the central findings (the multi-scale paradox, EMT-first
coordination, pathway selectivity, XBP1 consolidation, cross-platform pathway
concordance) are all confirmed by the data. The corrections are to specific numbers,
wordings, and one figure panel.

Status legend: **[done in repo]** applied to code/figures/tables here; **[author
action]** requires an edit to the manuscript or to a figure whose source is not in
this repository.

---

## A. Corrections required in the manuscript text

### A1. §3.2 — entropy values are medians, not a mean · [author action]

> Current: "CSC_TBXT+ cells showed the highest global entropy among all ten cell
> types, with **a mean of 9.97 bits**."

Every entropy value quoted in §3.2 matches the **median**, not the mean:

| Metric | Mean | Median | Quoted |
|---|---|---|---|
| global | 9.838 | 9.968 | 9.97 |
| ribosomal-excluded | 9.870 | 10.031 | 10.03 |
| immune evasion (norm) | 0.402 | 0.409 | 0.409 |
| ferroptosis (norm) | 0.433 | 0.438 | 0.438 |
| EMT (norm) | 0.414 | 0.437 | 0.437 |
| housekeeping (norm) | 0.683 | 0.683 | 0.683 |

Figure 2A already reads "ranked by median entropy." Change "a mean of 9.97 bits" to
"a median of 9.97 bits." The conclusion is unaffected: CSC ranks highest on either
statistic. *(Guarded by `test_supplementary_tables.py::TestS3EntropyStatistics`.)*

### A2. §2.4 — immune evasion is 28 genes, not 29 · [author action]

> Current: "immune evasion (**29 genes**)"

The set listed `PVRL2` and `NECTIN2`, which are two symbols for the same gene (PVRL2
is the withdrawn HGNC alias of NECTIN2). It only ever resolved to 28 distinct genes,
and every analysis output already reports 28. `PVRL2` has been removed from
`msep.pathways.IMMUNE_EVASION` **[done in repo]**; the text should read "28 genes,"
and Supplementary Table S1 now lists 28. *(Guarded by
`test_curation_and_cv.py::TestImmuneEvasionGeneSet`.)*

### A3. §3.4 — EMT ranking order · [author action]

> Current: "followed by basal cell carcinoma (5.49), osteosarcoma (5.25), and breast
> cancer (5.49)."

Listed out of order. Ascending EMT CV among the 10x cohorts after chordoma (4.632)
is osteosarcoma **5.25**, then basal cell carcinoma **5.49**, then breast cancer
**5.49**. Reorder to "followed by osteosarcoma (5.25), basal cell carcinoma (5.49),
and breast cancer (5.49)." *(Guarded by
`test_pan_cancer_table.py::TestEmtCoordination`.)*

### A4. §3.7 — SLC7A11 is not conserved in scRNA-seq · [author action]

> Current: "Gene-level analysis confirmed that TBXT (0.620), SMAD3 (0.540), SMAD2
> (0.700), SLC7A11 (0.563), and VIM (0.699) were among the most conserved genes
> across six tumors, **consistent with their low CV in the scRNA-seq data**."

The five bulk CVs are correct, but the "consistent with … scRNA-seq" clause fails for
most of them. scRNA-seq CV and rank among the 23 genes:

| Gene | scRNA-seq CV | rank |
|---|---|---|
| TBXT | 1.02 | 1 / 23 |
| VIM | 1.74 | 10 / 23 |
| SMAD2 | 1.83 | 12 / 23 |
| SMAD3 | 2.05 | 14 / 23 |
| SLC7A11 | **13.73** | **23 / 23** |

SLC7A11 is the *most* variable gene across CSC, not among the most conserved. Only
TBXT genuinely supports the original sentence. Suggested replacement (also handles
A5 and the S15 pointer):

> Gene-level analysis identified TBXT (CV = 0.620), SMAD3 (0.540), SMAD2 (0.700),
> SLC7A11 (0.563), and VIM (0.699) as the most conserved genes across the six
> tumors. For TBXT this mirrors its rank as the least variable gene in the scRNA-seq
> data, whereas SLC7A11 is the most variable gene across CSC (CV = 13.73) despite its
> uniformity across tumors. Consistent with this, gene-level CV did not correlate
> across the two platforms (Spearman ρ = −0.110, p = 0.618, n = 23; Supplementary
> Figure S15). Cross-platform concordance in this dataset is therefore a property of
> aggregate pathway behaviour rather than of individual genes, as expected given that
> the two measurements quantify different sources of variation — cell-to-cell within
> one tumor versus tumor-to-tumor across patients.

*(Guarded by `test_figure6_concordance.py::TestSupplementaryGeneLevel`.)*

### A5. §3.7 — notochord EMT wording · [author action]

> Current: "notochord EMT CV (0.282) was the lowest among **all** notochord
> pathways."

Housekeeping is lower (0.207). EMT is the lowest among the three **resistance**
pathways, which is the comparison the argument needs and matches the resistance-only
framing §3.3 already uses for chordoma. Change "all notochord pathways" to "the
notochord resistance pathways." *(Guarded by
`test_figure6_concordance.py::TestPanelB`.)*

---

## B. Figure corrections

### B1. Figure 6 — panel D removed, gene-level scatter moved to S15 · [done in repo]

The published panel D reported `Spearman ρ = 0.000, p = 1.0000`. That statistic
matched neither the source data (23 genes, ρ = −0.110, p = 0.618) nor the eight
points the panel plotted (ρ = −0.500, p = 0.207), and several coordinates were wrong
(SLC7A11 at ~2.1 vs a true 13.73). Recomputed honestly the correlation is weak and
not significant, so it cannot carry a main-figure concordance claim.

Figure 6 is now three panels (A–C, `figures/Publication/Figure_6_BulkValidation`).
The gene-level scatter, recomputed over all 23 genes, is Supplementary Figure S15
(`figures/Supplementary/Figure_S15_GeneLevel_CV_Concordance`). The superseded
`Figure_6D_GeneLevel_CV_Concordance_REAL.*` can be deleted. Update the Figure 6
legend to three panels and add the S15 pointer (see A4 text). *(Rebuilt by
`notebooks/paper_figures/figure6_bulk_validation.py`.)*

### B2. Ferroptosis gene set standardized; text and Figure 3 updated · [Figure 3 = author action]

Resolved by re-running from the chordoma `.h5ad`. The discrepancy was not a 33-vs-32
gene issue but **two different 33-gene ferroptosis sets** sharing only 28 genes:

| | Genes |
|---|---|
| Analysis set (Figure 2, scVI, sensitivity) | includes HSPB1, IREB2, NFS1, PROM2, STEAP3 → CV **5.61** |
| Old package/pan-cancer set | includes ALOX15B, CARS1, CISD2, PEBP1, SLC3A2 → CV **5.387** |

All 33 analysis-set genes are present and expressed in CSC (confirmed on the h5ad), so
the **analysis set is canonical**. The package, the pan-cancer table, and Table S1 are
now standardized to it **[done in repo]**:

- `msep.pathways.FERROPTOSIS` swapped to the analysis set (5 genes).
- `data/pan_cancer_cv.csv`: chordoma ferroptosis 5.387 → **5.609**; rank **5th → 7th**
  of twelve (robust: the types straddling chordoma are gene-set-invariant).
- Table S1 regenerated with the correct genes.
- Text edits applied to `manuscript/main.tex`: §3.3 "second-lowest CV (**5.61**)",
  §3.4 "**seventh** of twelve (CV = **5.61**)".

The 7th-place rank is still moderate coordination, well inside the solid-tumor range;
the pathway-selectivity thesis (EMT far more coordinated than ferroptosis or immune
evasion) is unaffected. *(Guarded by `test_curation_and_cv.py::TestFerroptosisGeneSet`,
`::TestFerroptosisCanonicalValue`, and `test_pan_cancer_table.py`.)*

**Figure 2** shows 5.61 and is therefore already correct — no change.

**Figure 3 · [author action]:** the published `Figure_3_PanCancer.png` shows chordoma
ferroptosis **5.39** and must be regenerated to **5.61** (Panel A heatmap cell + Panel C
ferroptosis bar). The figure source is in Colab; regenerate the chordoma column with the
canonical set. For full visual consistency, regenerating all 12 cancer types' ferroptosis
column on the canonical set is ideal, but only chordoma is quoted in the text and the
rank is settled. (The tracked reproducibility notebook `figure3_pan_cancer.ipynb` already
renders 5.61 from `data/pan_cancer_cv.csv`, in a different layout from the published
figure.)

---

## C. Supplementary deliverables

### C1. Tables S1–S10 assembled · [done in repo]

None existed as files before. All ten are now built by
`notebooks/paper_figures/build_supplementary_tables.py` into
`figures/Supplementary/tables/` (one CSV each plus a combined workbook). Two columns
need author input and are emitted empty: **S1 `literature_source`** (per-gene
citations; Methods §2.4 records only pathway-level provenance) and **S2 `reference`**.

### C2. Figure exports still needed · [author action]

- **S2** (ribosomal correction impact) and **S3** (entropy/CV UMAPs) exist as content
  under `figures/03_Entropy/` but are not exported under their S-numbers. S3 as
  declared also promises an across-cells-CV UMAP panel, which cannot be drawn per
  cell (CV is a population quantity); either drop that clause from the legend or add
  a per-cell-type CV panel.
- **S5** and **S6** currently share one file (`Figure_S5_S6_PerPatient_NK.pdf`); most
  journals want them split.
- Missing PNG companions: S1a, S1b, S1c, S5/S6, S7.

---

## Guard-test index

| Claim | Test |
|---|---|
| Pan-cancer ranks & EMT order (A3) | `test_pan_cancer_table.py` |
| Figure 6 A–C + gene-level null (A4, A5, B1) | `test_figure6_concordance.py` |
| Supplementary table shapes & §3.2 medians (A1) | `test_supplementary_tables.py` |
| PVRL2/NECTIN2, gene count, ferroptosis 5.387 (A2, B2) | `test_curation_and_cv.py` |

Run `pytest` to confirm the current tree still satisfies every pinned claim (122
tests at the time of writing).
