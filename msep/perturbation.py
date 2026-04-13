"""
In silico pseudo-perturbation analysis.

Natural-variance pseudo-perturbation stratifies cells by expression of
a *shield gene* (top vs. bottom quartile) and tests whether high
expression of one gene alters the population CV of other pathways.

XBP1 stress consolidation analysis stratifies cells into XBP1-zero,
XBP1-low, and XBP1-high groups and compares pathway coordination.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix

from .coordination import pathway_cv

__all__ = [
    "pseudo_perturbation",
    "xbp1_consolidation",
    "stratify_by_gene",
]


def stratify_by_gene(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    gene: str,
    lower_q: float = 0.25,
    upper_q: float = 0.75,
) -> dict:
    """Split cells into high/low groups based on a single gene.

    Only cells with nonzero expression are considered.  The thresholds
    are computed on the nonzero distribution.

    Parameters
    ----------
    count_matrix
        Cells × genes (raw counts).
    var_names
        Gene names.
    gene
        Gene symbol to stratify by.
    lower_q, upper_q
        Quantile thresholds on the nonzero expression distribution.

    Returns
    -------
    dict
        Keys ``"high_mask"``, ``"low_mask"`` (boolean arrays, shape
        n_cells), ``"n_high"``, ``"n_low"``, ``"threshold_high"``,
        ``"threshold_low"``.
    """
    idx = _gene_index(var_names, gene)
    col = count_matrix[:, idx]
    if issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()

    nonzero_mask = col > 0
    if nonzero_mask.sum() < 10:
        return {
            "high_mask": np.zeros(len(col), dtype=bool),
            "low_mask": np.zeros(len(col), dtype=bool),
            "n_high": 0, "n_low": 0,
            "threshold_high": np.nan, "threshold_low": np.nan,
        }

    nonzero_vals = col[nonzero_mask]
    thresh_low = np.percentile(nonzero_vals, 100 * lower_q)
    thresh_high = np.percentile(nonzero_vals, 100 * upper_q)

    high_mask = nonzero_mask & (col >= thresh_high)
    low_mask = nonzero_mask & (col <= thresh_low)

    return {
        "high_mask": high_mask,
        "low_mask": low_mask,
        "n_high": int(high_mask.sum()),
        "n_low": int(low_mask.sum()),
        "threshold_high": float(thresh_high),
        "threshold_low": float(thresh_low),
    }


def pseudo_perturbation(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    shield_genes: List[str],
    pathways: Dict[str, List[str]],
    n_perm: int = 500,
    seed: int = 42,
    lower_q: float = 0.25,
    upper_q: float = 0.75,
    min_genes: int = 3,
) -> pd.DataFrame:
    """Test cross-pathway coordination via natural-variance stratification.

    For each *shield_gene*, cells are split into high (top quartile) and
    low (bottom quartile) expressors.  Pathway CV is computed for each
    group, and the difference (Δ CV = high − low) is tested by permutation.

    Parameters
    ----------
    count_matrix
        Cells × genes (raw counts).
    var_names
        Gene names.
    shield_genes
        Genes to stratify cells by (one at a time).
    pathways
        ``{pathway_name: [gene, …]}``.
    n_perm
        Permutation iterations for significance testing.
    seed
        Random seed.
    lower_q, upper_q
        Quartile thresholds passed to :func:`stratify_by_gene`.
    min_genes
        Forwarded to :func:`pathway_cv`.

    Returns
    -------
    DataFrame
        Columns: ``shield_gene``, ``target_pathway``, ``cv_high``,
        ``cv_low``, ``delta_cv``, ``perm_p``, ``n_high``, ``n_low``.
    """
    rng = np.random.RandomState(seed)
    rows = []

    for sg in shield_genes:
        strat = stratify_by_gene(
            count_matrix, var_names, sg, lower_q, upper_q
        )
        if strat["n_high"] < 10 or strat["n_low"] < 10:
            continue

        mat_high = count_matrix[strat["high_mask"]]
        mat_low = count_matrix[strat["low_mask"]]

        for pw_name, pw_genes in pathways.items():
            cv_high, _ = pathway_cv(mat_high, var_names, pw_genes, min_genes)
            cv_low, _ = pathway_cv(mat_low, var_names, pw_genes, min_genes)

            if np.isnan(cv_high) or np.isnan(cv_low):
                continue

            delta = cv_high - cv_low

            # Permutation test
            combined_mask = strat["high_mask"] | strat["low_mask"]
            combined_idx = np.where(combined_mask)[0]
            n_high = strat["n_high"]
            perm_deltas = []

            for _ in range(n_perm):
                shuffled = rng.permutation(combined_idx)
                perm_high = np.zeros(count_matrix.shape[0], dtype=bool)
                perm_low = np.zeros(count_matrix.shape[0], dtype=bool)
                perm_high[shuffled[:n_high]] = True
                perm_low[shuffled[n_high:]] = True

                cv_ph, _ = pathway_cv(
                    count_matrix[perm_high], var_names, pw_genes, min_genes
                )
                cv_pl, _ = pathway_cv(
                    count_matrix[perm_low], var_names, pw_genes, min_genes
                )
                if not (np.isnan(cv_ph) or np.isnan(cv_pl)):
                    perm_deltas.append(cv_ph - cv_pl)

            if perm_deltas:
                perm_deltas = np.array(perm_deltas)
                p_val = float(np.mean(np.abs(perm_deltas) >= np.abs(delta)))
            else:
                p_val = np.nan

            rows.append({
                "shield_gene": sg,
                "target_pathway": pw_name,
                "cv_high": round(cv_high, 4),
                "cv_low": round(cv_low, 4),
                "delta_cv": round(delta, 4),
                "perm_p": round(p_val, 6),
                "n_high": strat["n_high"],
                "n_low": strat["n_low"],
            })

    return pd.DataFrame(rows)


def xbp1_consolidation(
    count_matrix: Union[np.ndarray, spmatrix],
    var_names,
    pathways: Dict[str, List[str]],
    high_pct: float = 85.0,
    min_genes: int = 3,
) -> pd.DataFrame:
    """XBP1 stress-induced defense consolidation analysis.

    Cells are grouped into XBP1-zero, XBP1-low (nonzero, below
    *high_pct* percentile), and XBP1-high (≥ *high_pct* percentile of
    nonzero values).  Pathway CV is compared across groups.

    Parameters
    ----------
    count_matrix
        Cells × genes (raw counts).
    var_names
        Gene names.
    pathways
        ``{pathway_name: [gene, …]}``.
    high_pct
        Percentile threshold on nonzero XBP1 expression for the
        "high" group.
    min_genes
        Forwarded to :func:`pathway_cv`.

    Returns
    -------
    DataFrame
        Columns: ``xbp1_group``, ``pathway``, ``cv``, ``n_cells``.
    """
    idx = _gene_index(var_names, "XBP1")
    col = count_matrix[:, idx]
    if issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()

    zero_mask = col == 0
    nonzero_mask = col > 0
    nonzero_vals = col[nonzero_mask]

    if nonzero_mask.sum() < 10:
        return pd.DataFrame(
            columns=["xbp1_group", "pathway", "cv", "n_cells"]
        )

    thresh = np.percentile(nonzero_vals, high_pct)
    high_mask = nonzero_mask & (col >= thresh)
    low_mask = nonzero_mask & (col < thresh)

    groups = {
        "XBP1-zero": zero_mask,
        "XBP1-low": low_mask,
        "XBP1-high": high_mask,
    }

    rows = []
    for grp_name, grp_mask in groups.items():
        n = int(grp_mask.sum())
        if n < 5:
            continue
        mat_grp = count_matrix[grp_mask]
        for pw_name, pw_genes in pathways.items():
            cv, _ = pathway_cv(mat_grp, var_names, pw_genes, min_genes)
            rows.append({
                "xbp1_group": grp_name,
                "pathway": pw_name,
                "cv": round(cv, 4) if not np.isnan(cv) else np.nan,
                "n_cells": n,
            })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# internal
# -----------------------------------------------------------------------

def _gene_index(var_names, gene: str) -> int:
    import pandas as pd
    if not isinstance(var_names, pd.Index):
        var_names = pd.Index(var_names)
    locs = np.where(var_names == gene)[0]
    if len(locs) == 0:
        raise KeyError(f"Gene {gene!r} not found in var_names")
    return int(locs[0])
