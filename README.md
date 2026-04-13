# msep — Multi-Scale Entropy Profiling

[![PyPI version](https://img.shields.io/pypi/v/msep?color=blue)](https://pypi.org/project/msep/)
[![Python versions](https://img.shields.io/pypi/pyversions/msep?color=blue)](https://pypi.org/project/msep/)
[![Tests](https://img.shields.io/badge/tests-27%2F27-brightgreen)](tests/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Multi-Scale Entropy Profiling** integrates per-cell Shannon entropy with across-cells coefficient of variation (CV), decomposed by biological pathway, to characterise population-level transcriptomic coordination in single-cell RNA-seq data.

The framework reveals states invisible to single-scale analyses — such as populations that are *individually diverse* (high per-cell entropy) yet *collectively disciplined* (low across-cells CV) — and is applicable to any cancer type or cellular system.

> Çavuş & Kuşkucu (2026). *Multi-Scale Entropy Profiling Reveals Pathway-Selective Defense Coordination Across Cancer Types.*

---

## Installation

```bash
pip install msep
```

For plotting support:
```bash
pip install msep[plotting]
```

## Quick Start

```python
import scanpy as sc
import msep

# Load your annotated AnnData (must have raw counts)
adata = sc.read_h5ad("my_data.h5ad")

# Run full multi-scale entropy profiling
result = msep.profile(
    adata,
    pathways="cancer_defense",       # built-in: ferroptosis, immune evasion, EMT, housekeeping
    cell_type_key="cell_type",       # column in adata.obs
    layer="raw_counts",              # layer with unnormalized counts
)

# Inspect the multi-scale paradox
result.paradox_summary               # per-cell entropy vs population CV per cell type

# Pathway-level coordination table
result.pathway_cv                    # DataFrame: cell_type × pathway × CV

# Publication-ready figures
fig = msep.plot_paradox(result, cv_pathway="emt")
fig.savefig("paradox.pdf")
```

## Full Pipeline (including perturbation analysis)

```python
result = msep.profile(
    adata,
    pathways="cancer_defense",
    cell_type_key="cell_type",
    layer="raw_counts",
    # Enable all analyses
    compute_bootstrap=True,          # 95% CI for CV
    compute_gene_cv=True,            # per-gene CV
    compute_perturbation=True,       # cross-pathway coordination
    compute_xbp1=True,               # stress consolidation
    shield_genes=["VIM", "GPX4", "FTH1", "HLA-E", "B2M", "TBXT"],
    perturbation_cell_type="CSC",    # restrict to specific cell type
    n_boot=1000,
    n_perm=500,
)

# Cross-pathway shield coordination
result.perturbation                  # shield gene → target pathway Δ CV + p-value

# XBP1 stress consolidation
result.xbp1                          # CV by XBP1 group (zero/low/high)
result.consolidation_score()         # how many pathways consolidate

# Bootstrap confidence intervals
result.bootstrap                     # dict of CI results per cell_type|pathway

# Per-gene coordination
result.gene_cv["emt"]                # gene-level CV for EMT pathway
```

## Custom Pathways

```python
my_pathways = {
    "glycolysis": ["HK1", "HK2", "PKM", "LDHA", "ENO1"],
    "apoptosis": ["BCL2", "BAX", "CASP3", "CASP9", "TP53"],
    "stemness": ["SOX2", "POU5F1", "NANOG", "KLF4", "MYC"],
}

result = msep.profile(adata, pathways=my_pathways, cell_type_key="cell_type")
```

## What It Computes

| Scale | Metric | Question answered |
|-------|--------|-------------------|
| Per-cell | Shannon entropy | How many gene programs does this cell engage? |
| Per-cell × pathway | Pathway entropy | Is this cell focused or broad within each pathway? |
| Across-cells | Pathway CV | Do all cells express this pathway at the same level? |
| Across-cells | Bootstrap CI | Is the CV difference statistically robust? |
| Across-cells | Fano factor | Is CV driven by mean expression level? |
| Cross-pathway | Pseudo-perturbation | Does high expression of gene X tighten pathway Y? |
| Stress response | XBP1 consolidation | Does ER stress coordinate all defense programs? |

## Key Concepts

**Individually diverse, collectively disciplined:** A population where each cell engages many gene programs (high per-cell entropy) but all cells converge on the same expression levels (low across-cells CV). This cannot be detected by single-scale entropy analysis.

**Pathway-selective coordination:** Not all pathways are equally coordinated. EMT genes may be tightly locked while immune evasion genes are heterogeneous — revealing distinct defense architectures.

**Stress-induced consolidation:** Under XBP1-mediated ER stress, defense pathways can become simultaneously more coordinated — a phenomenon observed in multiple cancer types.

## Bayesian Validation (optional)

Validate CV-based findings using scVI's generative model, which
separates technical variance (dropout, library size) from biological
variance:

```bash
pip install msep[bayesian]   # adds scvi-tools + torch
```

```python
bayes = msep.bayesian_validate(
    adata,
    pathways=msep.pathways.CANCER_DEFENSE,
    cell_type_key="cell_type",
    cell_types=["CSC_TBXT+", "T_cell", "Macrophage"],
    batch_key="patient_id",
    n_posterior_samples=25,
)

bayes.table                  # raw CV, denoised CV, BDR per cell_type × pathway
bayes.concordance            # ranking concordance per cell type
bayes.is_concordant          # True if all rankings preserved after denoising
bayes.csc_summary("CSC")     # detailed view for a specific cell type
```

The Bayesian Dispersion Ratio (BDR = biological variance / total variance)
quantifies what fraction of observed variance is genuine biology vs
technical noise, providing model-based evidence for coordination claims.

## API Reference

### Core

- `msep.profile(adata, ...)` → `MSEPResult` — main entry point

### Low-level functions

- `msep.entropy.per_cell_entropy(count_matrix, ...)` → entropy, n_expressed
- `msep.entropy.normalized_entropy(entropy, n_expressed)` → normalized
- `msep.coordination.pathway_cv(matrix, var_names, genes)` → cv, n_genes
- `msep.coordination.pathway_cv_table(matrix, var_names, pathways, labels)` → DataFrame
- `msep.coordination.bootstrap_cv(matrix, var_names, genes, n_boot)` → dict
- `msep.coordination.fano_factor(matrix, var_names, genes)` → fano, n_genes
- `msep.coordination.gene_level_cv(matrix, var_names, genes)` → DataFrame
- `msep.perturbation.pseudo_perturbation(matrix, var_names, shield_genes, pathways)` → DataFrame
- `msep.perturbation.xbp1_consolidation(matrix, var_names, pathways)` → DataFrame

### Bayesian (requires scvi-tools)

- `msep.bayesian_validate(adata, pathways, ...)` → `BayesianResult`
- `BayesianResult.table` — full results DataFrame
- `BayesianResult.concordance` — ranking concordance per cell type
- `BayesianResult.is_concordant` — bool: all rankings preserved?
- `BayesianResult.csc_summary(cell_type)` — filtered view

### Plotting

- `msep.plot_entropy_violin(result)` → Figure
- `msep.plot_pathway_cv_heatmap(result)` → Figure
- `msep.plot_paradox(result, cv_pathway)` → Figure
- `msep.plot_pan_cancer(cv_data, pathway, highlight)` → Figure

## Citation

If you use msep in your research, please cite:

```bibtex
@article{cavus2026msep,
  title={Multi-Scale Entropy Profiling Reveals Pathway-Selective Defense 
         Coordination Across Cancer Types},
  author={{\c{C}}avu{\c{s}}, {\"O}zge A. and Ku{\c{s}}kucu, Ay{\c{s}}eg{\"u}l},
  year={2026},
  journal={[submitted]},
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
