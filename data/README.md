# Tracked analysis tables

Small, derived tables that figure notebooks read directly, so a reviewer can render a
figure without re-running the pipeline that produced the numbers.

These are **distillations, not raw data**. The pipeline output they come from lives in
`processed/`, which is large and untracked.

## `pan_cancer_cv.csv`

Across-cells coefficient of variation per pathway for the malignant compartment of each
cancer type in the pan-cancer panel — the source of **Figure 3** and Supplementary Table S5.

| Column | Meaning |
|---|---|
| `cancer_type` | Key used by the pipeline; stable across rebuilds |
| `display_name` | Label rendered in figures |
| `source` | `CellxGene Census`, `GEO`, or `This study` |
| `accession` | Census release, GEO accession, or repository |
| `platform` | `10x` or `Smart-seq2` |
| `n_cells` | Malignant cells profiled after subsampling |
| `cell_selection` | How the malignant compartment was defined |
| `in_primary_comparison` | `True` for 10x cohorts. Smart-seq2 melanoma is reported separately because of known platform effects on expression variance (Methods §2.1.3) |
| `<pathway>_cv` | Mean across-cells CV over detected pathway genes |
| `<pathway>_n_genes` | Pathway genes detected in that cohort |

### Provenance

```
examples/chordoma_msep.ipynb  (pan-cancer section)
    │   Census queries + GEO downloads + malignant filter + msep.profile
    ▼
processed/pan_cancer_census_cv_expanded.csv   (untracked, hours to produce)
    │   notebooks/pan_cancer/build_pan_cancer_table.py
    ▼
data/pan_cancer_cv.csv                        (tracked, seconds to produce)
    │
    ▼
notebooks/figures/figure3_pan_cancer.ipynb
```

### Rebuilding

```bash
python notebooks/pan_cancer/build_pan_cancer_table.py
pytest tests/test_pan_cancer_table.py
```

`tests/test_pan_cancer_table.py` asserts this table still supports every numeric claim in
Results §3.4 — the chordoma EMT minimum, its ferroptosis rank of 5/12, its immune evasion
rank of 9/12, and the quoted CV values for the neighbouring cancer types. A failure there
means the table and the manuscript have diverged and one of them needs updating.

Rank-based assertions depend on the panel size. When the panel is extended (WP-3 grows it to
20 cancer types), update the test and the corresponding sentence in Results §3.4 in the same
commit.
