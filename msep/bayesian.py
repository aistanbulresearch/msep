"""
Bayesian variance decomposition via scVI generative model.

Uses scVI (Lopez et al., Nature Methods 2018) to fit a hierarchical
negative-binomial generative model that explicitly separates technical
variance (library size, dropout, capture efficiency) from biological
variance (true expression heterogeneity).

This module provides model-based validation that pathway-level CV
rankings reflect genuine biological coordination rather than technical
artefacts.

Requires ``scvi-tools >= 1.1`` (install via ``pip install msep[bayesian]``).

References
----------
Lopez, R. et al. (2018). Deep generative modeling for single-cell
transcriptomics. Nature Methods, 15(12), 1053–1058.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

__all__ = [
    "bayesian_variance_decomposition",
    "BayesianResult",
]


class BayesianResult:
    """Container for Bayesian variance decomposition results.

    Attributes
    ----------
    table : DataFrame
        Per (cell_type, pathway) results with columns: ``cell_type``,
        ``pathway``, ``n_cells``, ``n_genes``, ``raw_cv``,
        ``denoised_cv``, ``bio_variance``, ``tech_variance``,
        ``total_variance``, ``BDR``, ``cv_change_pct``.
    concordance : DataFrame
        Per-cell-type ranking concordance (Spearman ρ, concordant bool).
    model_params : dict
        scVI training parameters used.
    """

    def __init__(self, table: pd.DataFrame, concordance: pd.DataFrame,
                 model_params: dict):
        self.table = table
        self.concordance = concordance
        self.model_params = model_params

    @property
    def is_concordant(self) -> bool:
        """True if pathway ranking is preserved for all cell types."""
        return bool(self.concordance["concordant"].all())

    def csc_summary(self, cell_type: str = "CSC_TBXT+") -> pd.DataFrame:
        """Return results for a specific cell type, sorted by raw CV."""
        return (
            self.table[self.table["cell_type"] == cell_type]
            .sort_values("raw_cv")
            .reset_index(drop=True)
        )

    def __repr__(self) -> str:
        n_ct = self.table["cell_type"].nunique()
        n_pw = self.table["pathway"].nunique()
        conc = self.concordance["concordant"].sum()
        total = len(self.concordance)
        return (
            f"BayesianResult(cell_types={n_ct}, pathways={n_pw}, "
            f"concordant={conc}/{total})"
        )


def bayesian_variance_decomposition(
    adata,
    pathways: Dict[str, List[str]],
    cell_type_key: str = "cell_type",
    cell_types: Optional[List[str]] = None,
    batch_key: Optional[str] = None,
    layer: Optional[str] = None,
    # scVI parameters
    n_latent: int = 30,
    n_layers: int = 2,
    max_epochs: int = 200,
    early_stopping: bool = True,
    early_stopping_patience: int = 15,
    lr: float = 1e-3,
    # posterior sampling
    n_posterior_samples: int = 25,
    library_size: float = 1e4,
    # HVG selection
    n_top_genes: int = 3000,
    # misc
    seed: int = 42,
) -> BayesianResult:
    """Run Bayesian variance decomposition using scVI.

    Fits an scVI negative-binomial generative model, samples from the
    posterior to obtain denoised expression estimates, then computes
    pathway-level biological vs technical variance and compares the
    denoised CV ranking against the raw CV ranking.

    Parameters
    ----------
    adata
        AnnData with raw counts.
    pathways
        ``{pathway_name: [gene, …]}``.
    cell_type_key
        Column in ``adata.obs`` with cell-type labels.
    cell_types
        Subset of cell types to analyse.  If ``None``, uses all.
    batch_key
        Column for batch correction (e.g. ``"patient_id"``).
    layer
        Layer with raw counts.  Auto-detected if ``None``.
    n_latent, n_layers
        scVI architecture parameters.
    max_epochs, early_stopping, early_stopping_patience, lr
        Training parameters.
    n_posterior_samples
        Number of posterior draws for variance estimation.
    library_size
        Target library size for normalised expression.
    n_top_genes
        Number of HVGs for scVI (pathway genes always included).
    seed
        Random seed.

    Returns
    -------
    BayesianResult
        Results including raw/denoised CV, BDR, and concordance.

    Raises
    ------
    ImportError
        If ``scvi-tools`` is not installed.
    """
    try:
        import scvi as _scvi
    except ImportError:
        raise ImportError(
            "scvi-tools is required for Bayesian analysis. "
            "Install with: pip install msep[bayesian]"
        )

    import scanpy as sc
    from scipy.sparse import issparse, csr_matrix
    from scipy.stats import spearmanr

    np.random.seed(seed)

    # ------------------------------------------------------------------
    # 1. Prepare data
    # ------------------------------------------------------------------
    if cell_types is not None:
        mask = adata.obs[cell_type_key].isin(cell_types)
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata.copy()
        cell_types = list(adata_sub.obs[cell_type_key].unique())

    # Get raw counts
    if layer is not None and layer in adata_sub.layers:
        adata_sub.X = adata_sub.layers[layer].copy()
    elif "raw_counts" in adata_sub.layers:
        adata_sub.X = adata_sub.layers["raw_counts"].copy()

    # Ensure integer counts
    if issparse(adata_sub.X):
        adata_sub.X = adata_sub.X.toarray()
    adata_sub.X = np.round(adata_sub.X).astype(int)
    adata_sub.X = csr_matrix(adata_sub.X)

    # Collect all pathway genes present in the dataset
    all_pw_genes = set()
    for genes in pathways.values():
        for g in genes:
            if g in adata_sub.var_names:
                all_pw_genes.add(g)

    # HVG selection (force-include pathway genes)
    sc.pp.highly_variable_genes(
        adata_sub, n_top_genes=n_top_genes, flavor="seurat_v3",
        batch_key=batch_key, subset=False,
    )
    for g in all_pw_genes:
        if g in adata_sub.var_names:
            adata_sub.var.loc[g, "highly_variable"] = True

    adata_sub = adata_sub[:, adata_sub.var["highly_variable"]].copy()

    # ------------------------------------------------------------------
    # 2. Fit scVI
    # ------------------------------------------------------------------
    setup_kwargs = {}
    if batch_key is not None:
        setup_kwargs["batch_key"] = batch_key

    _scvi.model.SCVI.setup_anndata(adata_sub, **setup_kwargs)

    model = _scvi.model.SCVI(
        adata_sub,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood="nb",
        dispersion="gene",
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        plan_kwargs={"lr": lr},
        check_val_every_n_epoch=5,
    )

    # ------------------------------------------------------------------
    # 3. Posterior sampling
    # ------------------------------------------------------------------
    denoised = model.get_normalized_expression(
        library_size=library_size, return_numpy=True,
    )

    posterior_samples = []
    for _ in range(n_posterior_samples):
        sample = model.get_normalized_expression(
            library_size=library_size, return_numpy=True, n_samples=1,
        )
        posterior_samples.append(sample)

    posterior_stack = np.stack(posterior_samples, axis=0)
    posterior_var = posterior_stack.var(axis=0)

    # ------------------------------------------------------------------
    # 4. Compute pathway metrics
    # ------------------------------------------------------------------
    results = []
    ct_labels = np.asarray(adata_sub.obs[cell_type_key])

    for ct in cell_types:
        ct_mask = ct_labels == ct
        ct_idx = np.where(ct_mask)[0]
        n_cells = len(ct_idx)
        if n_cells < 5:
            continue

        for pw_name, pw_genes in pathways.items():
            pw_idx = []
            for g in pw_genes:
                if g in adata_sub.var_names:
                    pw_idx.append(list(adata_sub.var_names).index(g))

            if len(pw_idx) < 3:
                continue
            pw_idx = np.array(pw_idx)

            # Raw CV
            raw_X = adata_sub.X[ct_idx][:, pw_idx]
            if issparse(raw_X):
                raw_X = raw_X.toarray()
            raw_X = raw_X.astype(float)
            raw_means = raw_X.mean(axis=0)
            raw_stds = raw_X.std(axis=0)
            valid_raw = raw_means > 0
            raw_cv = (
                float(np.mean(raw_stds[valid_raw] / raw_means[valid_raw]))
                if valid_raw.sum() > 0 else np.nan
            )

            # Denoised CV
            den_X = denoised[ct_idx][:, pw_idx]
            den_means = den_X.mean(axis=0)
            den_stds = den_X.std(axis=0)
            valid_den = den_means > 0
            denoised_cv = (
                float(np.mean(den_stds[valid_den] / den_means[valid_den]))
                if valid_den.sum() > 0 else np.nan
            )

            # Biological variance
            bio_var = den_X.var(axis=0)
            bio_var_mean = float(np.mean(bio_var[valid_den]))

            # Technical variance
            tech_var = posterior_var[ct_idx][:, pw_idx].mean(axis=0)
            tech_var_mean = float(np.mean(tech_var[valid_den]))

            # BDR
            total_var = bio_var_mean + tech_var_mean
            bdr = bio_var_mean / total_var if total_var > 0 else np.nan

            cv_change = (
                round(100 * (denoised_cv - raw_cv) / raw_cv, 1)
                if raw_cv > 0 and not np.isnan(raw_cv) else np.nan
            )

            results.append({
                "cell_type": ct,
                "pathway": pw_name,
                "n_cells": n_cells,
                "n_genes": len(pw_idx),
                "raw_cv": round(raw_cv, 4),
                "denoised_cv": round(denoised_cv, 4),
                "bio_variance": round(bio_var_mean, 6),
                "tech_variance": round(tech_var_mean, 6),
                "total_variance": round(total_var, 6),
                "BDR": round(bdr, 4),
                "cv_change_pct": cv_change,
            })

    df_results = pd.DataFrame(results)

    # ------------------------------------------------------------------
    # 5. Concordance
    # ------------------------------------------------------------------
    concordance_rows = []
    for ct in cell_types:
        ct_data = df_results[df_results["cell_type"] == ct]
        if len(ct_data) < 2:
            continue
        raw_rank = ct_data.sort_values("raw_cv")["pathway"].tolist()
        den_rank = ct_data.sort_values("denoised_cv")["pathway"].tolist()
        rho, _ = spearmanr(ct_data["raw_cv"], ct_data["denoised_cv"])
        concordance_rows.append({
            "cell_type": ct,
            "concordant": raw_rank == den_rank,
            "spearman_rho": round(rho, 4),
            "raw_ranking": " < ".join(raw_rank),
            "denoised_ranking": " < ".join(den_rank),
        })

    df_concordance = pd.DataFrame(concordance_rows)

    model_params = {
        "n_latent": n_latent, "n_layers": n_layers,
        "gene_likelihood": "nb", "dispersion": "gene",
        "max_epochs": max_epochs, "lr": lr,
        "n_posterior_samples": n_posterior_samples,
        "library_size": library_size,
        "n_top_genes": n_top_genes,
        "epochs_trained": model.history["elbo_train"].shape[0],
    }

    return BayesianResult(df_results, df_concordance, model_params)
