"""
Across-cells population coordination metrics.

The coefficient of variation (CV = σ/μ) of each gene is computed across
all cells of a given type, then averaged over genes within a pathway.
Lower pathway CV indicates tighter population-level coordination.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix

__all__ = [
    "pathway_cv",
    "pathway_cv_table",
    "bootstrap_cv",
    "fano_factor",
    "gene_level_cv",
]


# -----------------------------------------------------------------------
# Core CV
# -----------------------------------------------------------------------

def pathway_cv(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    gene_list: List[str],
    min_genes: int = 3,
) -> tuple[float, int]:
    """Compute mean across-cells CV for a single pathway.

    Parameters
    ----------
    count_matrix
        Cells × genes (raw counts recommended; log-normalized accepted).
    var_names
        Gene names matching the columns of *count_matrix*.
    gene_list
        Genes belonging to the pathway.
    min_genes
        Minimum number of detected pathway genes required.

    Returns
    -------
    cv : float
        Mean CV across expressed pathway genes (``NaN`` if too few genes).
    n_genes : int
        Number of pathway genes found and used.
    """
    present, _, indices = _resolve(var_names, gene_list)
    if len(present) < min_genes:
        return np.nan, len(present)

    sub = count_matrix[:, indices]
    if issparse(sub):
        sub = sub.toarray()
    sub = sub.astype(np.float64, copy=False)

    means = sub.mean(axis=0)
    stds = sub.std(axis=0)
    valid = means > 0
    if valid.sum() == 0:
        return np.nan, len(present)

    cvs = stds[valid] / means[valid]
    return float(np.mean(cvs)), len(present)


# -----------------------------------------------------------------------
# Tabular helpers
# -----------------------------------------------------------------------

def pathway_cv_table(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    pathways: Dict[str, List[str]],
    cell_type_labels: Optional[np.ndarray] = None,
    min_genes: int = 3,
) -> pd.DataFrame:
    """Compute pathway CV for every (cell-type, pathway) combination.

    Parameters
    ----------
    count_matrix
        Cells × genes.
    var_names
        Gene names.
    pathways
        ``{pathway_name: [gene, …]}``.
    cell_type_labels
        Per-cell labels.  If ``None``, all cells are treated as one group.
    min_genes
        Forwarded to :func:`pathway_cv`.

    Returns
    -------
    DataFrame
        Columns: ``cell_type``, ``pathway``, ``cv``, ``n_genes``.
    """
    rows = []
    if cell_type_labels is None:
        cell_type_labels = np.full(count_matrix.shape[0], "all")

    for ct in np.unique(cell_type_labels):
        mask = cell_type_labels == ct
        mat_ct = count_matrix[mask]
        for pw_name, pw_genes in pathways.items():
            cv, ng = pathway_cv(mat_ct, var_names, pw_genes, min_genes)
            rows.append(
                {"cell_type": ct, "pathway": pw_name, "cv": cv, "n_genes": ng}
            )
    return pd.DataFrame(rows)


def gene_level_cv(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    gene_list: List[str],
) -> pd.DataFrame:
    """Per-gene CV across cells.

    Returns
    -------
    DataFrame
        Columns: ``gene``, ``mean``, ``std``, ``cv``.
    """
    present, _, indices = _resolve(var_names, gene_list)
    sub = count_matrix[:, indices]
    if issparse(sub):
        sub = sub.toarray()
    sub = sub.astype(np.float64, copy=False)

    means = sub.mean(axis=0)
    stds = sub.std(axis=0)
    cvs = np.where(means > 0, stds / means, np.nan)

    return pd.DataFrame({
        "gene": present,
        "mean": means,
        "std": stds,
        "cv": cvs,
    })


# -----------------------------------------------------------------------
# Bootstrap CI
# -----------------------------------------------------------------------

def bootstrap_cv(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    gene_list: List[str],
    n_boot: int = 1000,
    ci: float = 0.95,
    seed: int = 42,
    min_genes: int = 3,
) -> dict:
    """Bootstrap confidence interval for pathway CV.

    Parameters
    ----------
    count_matrix
        Cells × genes.
    var_names
        Gene names.
    gene_list
        Pathway genes.
    n_boot
        Number of bootstrap iterations.
    ci
        Confidence level (default 0.95 → 2.5th–97.5th percentile).
    seed
        Random seed for reproducibility.
    min_genes
        Minimum pathway genes required.

    Returns
    -------
    dict
        Keys: ``cv_mean``, ``ci_low``, ``ci_high``, ``n_genes``,
        ``bootstrap_distribution``.
    """
    present, _, indices = _resolve(var_names, gene_list)
    if len(present) < min_genes:
        return {
            "cv_mean": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            "n_genes": len(present), "bootstrap_distribution": [],
        }

    sub = count_matrix[:, indices]
    if issparse(sub):
        sub = sub.toarray()
    sub = sub.astype(np.float64, copy=False)
    n_cells = sub.shape[0]

    rng = np.random.RandomState(seed)
    boot_cvs = []
    for _ in range(n_boot):
        idx = rng.choice(n_cells, n_cells, replace=True)
        sample = sub[idx]
        means = sample.mean(axis=0)
        stds = sample.std(axis=0)
        valid = means > 0
        if valid.sum() > 0:
            boot_cvs.append(float(np.mean(stds[valid] / means[valid])))

    alpha = (1 - ci) / 2
    return {
        "cv_mean": float(np.mean(boot_cvs)) if boot_cvs else np.nan,
        "ci_low": float(np.percentile(boot_cvs, 100 * alpha)) if boot_cvs else np.nan,
        "ci_high": float(np.percentile(boot_cvs, 100 * (1 - alpha))) if boot_cvs else np.nan,
        "n_genes": len(present),
        "bootstrap_distribution": boot_cvs,
    }


# -----------------------------------------------------------------------
# Fano factor (variance / mean) — mean-adjusted dispersion control
# -----------------------------------------------------------------------

def fano_factor(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    gene_list: List[str],
    min_genes: int = 3,
) -> tuple[float, int]:
    """Mean Fano factor (variance/mean) across pathway genes.

    Used as a control metric to confirm that CV differences are not
    driven solely by differences in mean expression level.

    Returns
    -------
    fano : float
        Mean Fano factor.
    n_genes : int
        Number of pathway genes used.
    """
    present, _, indices = _resolve(var_names, gene_list)
    if len(present) < min_genes:
        return np.nan, len(present)

    sub = count_matrix[:, indices]
    if issparse(sub):
        sub = sub.toarray()
    sub = sub.astype(np.float64, copy=False)

    means = sub.mean(axis=0)
    variances = sub.var(axis=0)
    valid = means > 0
    if valid.sum() == 0:
        return np.nan, len(present)

    fanos = variances[valid] / means[valid]
    return float(np.mean(fanos)), len(present)


# -----------------------------------------------------------------------
# internal
# -----------------------------------------------------------------------

def _resolve(var_names, gene_list):
    from .pathways import resolve_genes
    return resolve_genes(var_names, gene_list)
