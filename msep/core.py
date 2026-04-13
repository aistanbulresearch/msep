"""
High-level profiling API.

``msep.profile(adata, ...)`` is the main entry point.  It orchestrates
per-cell entropy, across-cells CV, bootstrap CI, pseudo-perturbation,
and XBP1 consolidation analysis in a single call.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from .entropy import per_cell_entropy, normalized_entropy
from .coordination import (
    pathway_cv,
    pathway_cv_table,
    bootstrap_cv,
    gene_level_cv,
)
from .perturbation import pseudo_perturbation, xbp1_consolidation
from .pathways import get_pathways, resolve_genes
from .result import MSEPResult

__all__ = ["profile"]


def profile(
    adata,
    pathways: Union[str, Dict[str, List[str]]] = "cancer_defense",
    cell_type_key: str = "cell_type",
    layer: Optional[str] = None,
    # what to compute
    compute_bootstrap: bool = True,
    compute_gene_cv: bool = True,
    compute_perturbation: bool = False,
    compute_xbp1: bool = False,
    # perturbation settings
    shield_genes: Optional[List[str]] = None,
    perturbation_cell_type: Optional[str] = None,
    n_perm: int = 500,
    # bootstrap settings
    n_boot: int = 1000,
    bootstrap_cell_types: Optional[List[str]] = None,
    # entropy settings
    exclude_ribosomal: bool = True,
    base: int = 2,
    seed: int = 42,
) -> MSEPResult:
    """Run multi-scale entropy profiling on an AnnData object.

    Parameters
    ----------
    adata
        Annotated data matrix (``anndata.AnnData``).  Must contain raw
        counts in either ``adata.X``, ``adata.raw.X``, or the layer
        specified by *layer*.
    pathways
        Either a built-in name (``"cancer_defense"``) or a custom dict
        ``{pathway_name: [gene_symbol, …]}``.
    cell_type_key
        Column in ``adata.obs`` that holds cell-type labels.
    layer
        Layer containing raw counts.  If ``None``, tries ``"raw_counts"``,
        then ``adata.raw.X``, then ``adata.X``.
    compute_bootstrap
        Compute bootstrap 95% CI for pathway CV.
    compute_gene_cv
        Compute per-gene CV for each pathway.
    compute_perturbation
        Run pseudo-perturbation analysis (requires *shield_genes*).
    compute_xbp1
        Run XBP1 stress consolidation analysis.
    shield_genes
        Genes to use in pseudo-perturbation (e.g. ``["VIM", "GPX4"]``).
    perturbation_cell_type
        Restrict perturbation analysis to this cell type.
    n_perm
        Permutation iterations for pseudo-perturbation.
    n_boot
        Bootstrap iterations for CV confidence intervals.
    bootstrap_cell_types
        Cell types to bootstrap (default: all).
    exclude_ribosomal
        Exclude RPL/RPS genes from global entropy.
    base
        Logarithm base for entropy (2 = bits).
    seed
        Random seed for reproducibility.

    Returns
    -------
    MSEPResult
        Container with all computed metrics.

    Examples
    --------
    >>> import msep
    >>> result = msep.profile(adata, cell_type_key="cell_type_fine")
    >>> result.paradox_summary
    >>> result.pathway_cv
    """
    # ------------------------------------------------------------------
    # 0. Resolve inputs
    # ------------------------------------------------------------------
    if isinstance(pathways, str):
        pw_dict = get_pathways(pathways)
    else:
        pw_dict = get_pathways(custom=pathways)

    count_mat = _get_raw_counts(adata, layer)
    var_names = adata.var_names
    ct_labels = np.asarray(adata.obs[cell_type_key])

    # Validate pathways against dataset
    pw_report = {}
    for pw_name, pw_genes in pw_dict.items():
        present, missing, _ = resolve_genes(var_names, pw_genes)
        pw_report[pw_name] = {
            "present": len(present), "total": len(pw_genes),
            "missing": missing,
        }

    # Ribosomal gene indices
    ribo_idx = None
    if exclude_ribosomal:
        ribo_mask = var_names.str.startswith(("RPS", "RPL", "MRPL", "MRPS"))
        ribo_idx = np.where(ribo_mask)[0] if ribo_mask.any() else None

    # ------------------------------------------------------------------
    # 1. Per-cell entropy
    # ------------------------------------------------------------------
    # Global
    h_global, n_expr_global = per_cell_entropy(
        count_mat, base=base, exclude_indices=ribo_idx
    )

    entropy_df = pd.DataFrame({
        "cell_type": ct_labels,
        "entropy_global": h_global,
        "n_expressed": n_expr_global,
    }, index=adata.obs_names)

    # Pathway-specific
    for pw_name, pw_genes in pw_dict.items():
        _, _, pw_idx = resolve_genes(var_names, pw_genes)
        if len(pw_idx) < 3:
            entropy_df[f"entropy_{pw_name}"] = np.nan
            entropy_df[f"entropy_{pw_name}_norm"] = np.nan
            continue
        h_pw, n_pw = per_cell_entropy(count_mat, gene_indices=pw_idx, base=base)
        h_pw_norm = normalized_entropy(h_pw, n_pw, base=base)
        entropy_df[f"entropy_{pw_name}"] = h_pw
        entropy_df[f"entropy_{pw_name}_norm"] = h_pw_norm

    # ------------------------------------------------------------------
    # 2. Across-cells CV
    # ------------------------------------------------------------------
    cv_df = pathway_cv_table(
        count_mat, var_names, pw_dict,
        cell_type_labels=ct_labels,
    )

    # ------------------------------------------------------------------
    # 3. Gene-level CV (optional)
    # ------------------------------------------------------------------
    gene_cv_dict = None
    if compute_gene_cv:
        gene_cv_dict = {}
        for pw_name, pw_genes in pw_dict.items():
            gene_cv_dict[pw_name] = gene_level_cv(
                count_mat, var_names, pw_genes
            )

    # ------------------------------------------------------------------
    # 4. Bootstrap CI (optional)
    # ------------------------------------------------------------------
    boot_dict = None
    if compute_bootstrap:
        boot_dict = {}
        cts = bootstrap_cell_types or list(np.unique(ct_labels))
        for ct in cts:
            mask = ct_labels == ct
            mat_ct = count_mat[mask]
            for pw_name, pw_genes in pw_dict.items():
                result = bootstrap_cv(
                    mat_ct, var_names, pw_genes,
                    n_boot=n_boot, seed=seed,
                )
                # drop full distribution to keep memory low
                result.pop("bootstrap_distribution", None)
                boot_dict[f"{ct}|{pw_name}"] = result

    # ------------------------------------------------------------------
    # 5. Pseudo-perturbation (optional)
    # ------------------------------------------------------------------
    perturb_df = None
    if compute_perturbation and shield_genes:
        if perturbation_cell_type is not None:
            mask = ct_labels == perturbation_cell_type
            mat_p = count_mat[mask]
        else:
            mat_p = count_mat
        perturb_df = pseudo_perturbation(
            mat_p, var_names, shield_genes, pw_dict,
            n_perm=n_perm, seed=seed,
        )

    # ------------------------------------------------------------------
    # 6. XBP1 consolidation (optional)
    # ------------------------------------------------------------------
    xbp1_df = None
    if compute_xbp1:
        if perturbation_cell_type is not None:
            mask = ct_labels == perturbation_cell_type
            mat_x = count_mat[mask]
        else:
            mat_x = count_mat
        try:
            xbp1_df = xbp1_consolidation(mat_x, var_names, pw_dict)
        except KeyError:
            xbp1_df = None  # XBP1 not in dataset

    # ------------------------------------------------------------------
    # 7. Assemble result
    # ------------------------------------------------------------------
    metadata = {
        "n_cells": count_mat.shape[0],
        "n_genes": count_mat.shape[1],
        "cell_types": list(np.unique(ct_labels)),
        "pathway_report": pw_report,
        "base": base,
        "exclude_ribosomal": exclude_ribosomal,
        "seed": seed,
    }

    return MSEPResult(
        per_cell_entropy=entropy_df,
        pathway_cv=cv_df,
        gene_cv=gene_cv_dict,
        bootstrap=boot_dict,
        perturbation=perturb_df,
        xbp1=xbp1_df,
        pathways_used=pw_dict,
        metadata=metadata,
    )


# -----------------------------------------------------------------------
# internal
# -----------------------------------------------------------------------

def _get_raw_counts(adata, layer: Optional[str]):
    """Extract raw count matrix from AnnData."""
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(
                f"Layer {layer!r} not found. Available: {list(adata.layers)}"
            )
        return adata.layers[layer]

    # auto-detect
    if "raw_counts" in adata.layers:
        return adata.layers["raw_counts"]
    if adata.raw is not None:
        return adata.raw.X
    return adata.X
