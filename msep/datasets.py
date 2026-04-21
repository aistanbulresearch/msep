"""
Built-in example datasets for quick demonstration of msep.

The main entry point is :func:`load_example`, which returns a small,
deterministic synthetic :class:`anndata.AnnData` object that reproduces the
qualitative chordoma MSEP paradox (high per-cell entropy in CSC combined with
low across-cells EMT CV and high immune-evasion CV) on a scale suitable for
a 5-minute Colab quickstart.

The dataset is generated on-the-fly from a seeded NumPy RNG; nothing is
bundled as a binary asset, so the package install remains small.

Examples
--------
>>> import msep
>>> adata = msep.datasets.load_example()
>>> adata
AnnData object with n_obs x n_vars = 500 x ...
    obs: 'cell_type', 'patient_id'
    layers: 'raw_counts'
>>> result = msep.profile(adata, pathways="cancer_defense")
>>> result.paradox_summary
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import anndata as ad

from .pathways import CANCER_DEFENSE

__all__ = ["load_example", "list_examples"]


# ---------------------------------------------------------------------------
# Per-(cell_type, pathway) expression specification
# ---------------------------------------------------------------------------

# (mean, sd_mean, dispersion)
#   mean       — base mean count across cells
#   sd_mean    — spread of the per-gene base mean WITHIN a (cell_type,pathway)
#                block (larger values => more gene-to-gene heterogeneity)
#   dispersion — NB-like over-dispersion sampled on top of the base mean
#                (larger values => more cell-to-cell heterogeneity per gene)
#
# These values are hand-tuned so that msep.profile() recovers the paper's
# qualitative ranking on the synthetic output:
#   CSC_TBXT+: EMT CV < ferroptosis CV < immune_evasion CV, and CSC has the
#   highest per-cell entropy of the three cell types.

_EXPR_SPEC: Dict[str, Dict[str, Tuple[float, float, float]]] = {
    "CSC_TBXT+": {
        "emt":            (8.0,  1.0, 0.25),
        "ferroptosis":    (5.0,  1.5, 0.55),
        "immune_evasion": (3.0,  2.5, 1.80),
        "housekeeping":   (12.0, 0.8, 0.18),
    },
    "T_cell": {
        "emt":            (1.5,  2.0, 1.40),
        "ferroptosis":    (2.0,  1.8, 1.10),
        "immune_evasion": (6.0,  2.5, 1.20),
        "housekeeping":   (10.0, 1.0, 0.22),
    },
    "Macrophage": {
        "emt":            (2.0,  2.0, 1.50),
        "ferroptosis":    (3.0,  1.8, 1.05),
        "immune_evasion": (5.0,  3.0, 1.60),
        "housekeeping":   (11.0, 0.9, 0.22),
    },
    "Stromal": {
        "emt":            (6.0,  2.0, 1.20),
        "ferroptosis":    (3.0,  1.8, 1.10),
        "immune_evasion": (2.0,  2.0, 1.60),
        "housekeeping":   (10.5, 1.0, 0.24),
    },
}

# Background ("non-pathway") genes — drawn from a broad exponential mean
# distribution to give every cell a rich transcriptomic background and so
# produce non-degenerate per-cell Shannon entropy.
_N_BACKGROUND_GENES = 120


def _sample_nb(rng: np.random.Generator, mean: np.ndarray,
               dispersion: float) -> np.ndarray:
    """
    Sample from a negative-binomial parametrised by (mean, dispersion).

    Variance = mean + dispersion * mean^2. When ``dispersion -> 0`` the
    distribution collapses to Poisson.
    """
    mean = np.asarray(mean, dtype=np.float64)
    mean = np.maximum(mean, 1e-6)
    if dispersion <= 0:
        return rng.poisson(mean).astype(np.float32)
    # NB parametrisation: n = 1/dispersion, p = n / (n + mean)
    n = 1.0 / dispersion
    p = n / (n + mean)
    return rng.negative_binomial(n, p).astype(np.float32)


# Extra single-gene markers that are not part of CANCER_DEFENSE but are
# required for downstream demos (XBP1 stress consolidation). Kept outside
# the pathway dict so pathway_cv remains uncontaminated.
_EXTRA_GENES: List[str] = ["XBP1"]


def _build_gene_universe() -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Union of the four cancer-defense pathways plus generic background genes
    plus a small set of extra markers (e.g. XBP1) required by downstream
    demos.

    Returns ``(gene_names, pathway_to_genes)``.
    """
    pathway_to_genes: Dict[str, List[str]] = {}
    seen: List[str] = []

    for pw_name, genes in CANCER_DEFENSE.items():
        in_pw: List[str] = []
        for g in genes:
            if g not in seen:
                seen.append(g)
            in_pw.append(g)
        pathway_to_genes[pw_name] = in_pw

    for g in _EXTRA_GENES:
        if g not in seen:
            seen.append(g)

    background = [f"BG_{i:03d}" for i in range(_N_BACKGROUND_GENES)]
    gene_names = seen + background
    return gene_names, pathway_to_genes


def _simulate_block(
    rng: np.random.Generator,
    n_cells: int,
    gene_base_means: np.ndarray,
    dispersion: float,
) -> np.ndarray:
    """
    Simulate ``n_cells x len(gene_base_means)`` counts for one
    (cell_type, pathway) block.
    """
    mean_matrix = np.broadcast_to(
        gene_base_means.reshape(1, -1),
        (n_cells, gene_base_means.size),
    )
    return _sample_nb(rng, mean_matrix, dispersion)


def load_example(
    n_cells: int = 500,
    seed: int = 42,
    n_patients: int = 3,
    include_xbp1_stress: bool = True,
) -> ad.AnnData:
    """Return a small synthetic :class:`anndata.AnnData` for quick demos.

    The output is deterministic given the seed and is designed to reproduce
    the paper's qualitative multi-scale paradox: the ``CSC_TBXT+`` cells
    exhibit the highest per-cell Shannon entropy among the included cell
    types, the lowest EMT coefficient of variation, and the highest
    immune-evasion coefficient of variation.

    Parameters
    ----------
    n_cells
        Total number of cells in the returned AnnData. Split across four
        cell types in fixed proportions (CSC 40%, T cell 30%, macrophage
        20%, stromal 10%). Must be >= 40.
    seed
        Seed for the internal NumPy RNG. Changing the seed produces a
        statistically equivalent but non-identical dataset.
    n_patients
        How many synthetic patients to assign cells to (uniform split per
        cell type). Useful for demonstrating the ``patient_id`` /
        ``batch_key`` machinery in ``msep.bayesian_validate``.
    include_xbp1_stress
        If True (default), the ``XBP1`` gene is set up so that roughly 45%
        of CSC cells have zero expression and the rest span a long-tailed
        non-zero distribution, making the dataset suitable for
        :func:`msep.perturbation.xbp1_consolidation` demos.

    Returns
    -------
    anndata.AnnData
        Shape ``(n_cells, ~235)`` with:

        * ``adata.X``                 — sparse CSR raw counts
        * ``adata.layers["raw_counts"]`` — duplicate for downstream tools
        * ``adata.obs["cell_type"]``   — categorical, four levels
        * ``adata.obs["patient_id"]``  — categorical, ``n_patients`` levels
        * ``adata.uns["pathways"]``    — the pathway dict used to simulate

    Notes
    -----
    This is a *toy* dataset intended for quickstart and unit tests only —
    do not draw biological conclusions from it. For the full chordoma
    single-cell pipeline see ``examples/chordoma_msep.ipynb``.
    """
    if n_cells < 40:
        raise ValueError(
            f"n_cells must be >= 40 to populate all four cell types "
            f"(got n_cells={n_cells})"
        )

    rng = np.random.default_rng(seed)

    # ---- cell-type composition ------------------------------------------------
    proportions = {
        "CSC_TBXT+": 0.40,
        "T_cell":    0.30,
        "Macrophage": 0.20,
        "Stromal":   0.10,
    }
    sizes: Dict[str, int] = {
        ct: int(round(n_cells * p)) for ct, p in proportions.items()
    }
    # correct rounding drift so the total exactly matches n_cells
    drift = n_cells - sum(sizes.values())
    sizes["CSC_TBXT+"] += drift

    # ---- gene universe --------------------------------------------------------
    gene_names, pathway_to_genes = _build_gene_universe()
    gene_index = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    # ---- base means for background genes (shared across cell types) ----------
    # Exponential-distributed background means give each cell a rich,
    # heterogeneous non-pathway transcriptome.
    bg_means = rng.exponential(scale=3.0, size=_N_BACKGROUND_GENES)
    bg_dispersion = 0.45

    # ---- simulate -------------------------------------------------------------
    blocks: List[np.ndarray] = []
    cell_type_labels: List[str] = []

    for cell_type, n_ct in sizes.items():
        spec = _EXPR_SPEC[cell_type]
        ct_matrix = np.zeros((n_ct, n_genes), dtype=np.float32)

        for pw_name, genes in pathway_to_genes.items():
            mean_base, mean_sd, dispersion = spec[pw_name]
            # per-gene base mean (shared across cells within this block)
            gene_base_means = np.maximum(
                rng.normal(loc=mean_base, scale=mean_sd, size=len(genes)),
                0.1,
            )
            block = _simulate_block(rng, n_ct, gene_base_means, dispersion)

            for j, g in enumerate(genes):
                ct_matrix[:, gene_index[g]] = block[:, j]

        # background genes
        bg_block = _simulate_block(rng, n_ct, bg_means, bg_dispersion)
        bg_slice = slice(n_genes - _N_BACKGROUND_GENES, n_genes)
        ct_matrix[:, bg_slice] = bg_block

        # extra single-gene markers (XBP1, …) get a generic low-baseline
        # expression here; cell-type-specific patterns are applied below.
        for extra in _EXTRA_GENES:
            if extra in gene_index:
                ct_matrix[:, gene_index[extra]] = rng.poisson(
                    1.5, size=n_ct
                ).astype(np.float32)

        blocks.append(ct_matrix)
        cell_type_labels.extend([cell_type] * n_ct)

    counts = np.vstack(blocks)

    # ---- TBXT addiction (low-CV within CSC) ----------------------------------
    if "TBXT" in gene_index:
        csc_mask = np.array(cell_type_labels) == "CSC_TBXT+"
        tbxt_idx = gene_index["TBXT"]
        # Non-CSC cells have essentially no TBXT.
        counts[~csc_mask, tbxt_idx] = rng.poisson(0.1, size=(~csc_mask).sum())
        # CSC cells: tightly-locked around a shared high mean (low CV).
        counts[csc_mask, tbxt_idx] = rng.normal(
            loc=18.0, scale=1.5, size=csc_mask.sum()
        ).clip(min=0).round()

    # ---- XBP1 stress pattern (nonzero/high only in a subset of CSC) ---------
    if include_xbp1_stress and "XBP1" in gene_index:
        csc_idx = np.where(np.array(cell_type_labels) == "CSC_TBXT+")[0]
        xbp1_idx = gene_index["XBP1"]
        # Default: most CSC zero
        counts[:, xbp1_idx] = 0
        # 43% non-zero (matches paper), of which the top ~9% are "high"
        rng.shuffle(csc_idx)
        n_nonzero = int(len(csc_idx) * 0.43)
        n_high = int(len(csc_idx) * 0.09)
        low_cells = csc_idx[:n_nonzero - n_high]
        high_cells = csc_idx[n_nonzero - n_high: n_nonzero]
        counts[low_cells, xbp1_idx] = rng.poisson(2.0, size=low_cells.size)
        counts[high_cells, xbp1_idx] = rng.poisson(15.0, size=high_cells.size)
        # Other cell types: normal background expression
        non_csc = np.where(np.array(cell_type_labels) != "CSC_TBXT+")[0]
        counts[non_csc, xbp1_idx] = rng.poisson(1.0, size=non_csc.size)

    # ---- assemble AnnData -----------------------------------------------------
    patient_ids = [f"P{(i % n_patients) + 1:02d}" for i in range(n_cells)]

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                cell_type_labels,
                categories=list(proportions.keys()),
            ),
            "patient_id": pd.Categorical(patient_ids),
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)

    sparse_counts = csr_matrix(counts)
    adata = ad.AnnData(X=sparse_counts, obs=obs, var=var)
    adata.layers["raw_counts"] = sparse_counts.copy()
    adata.uns["pathways"] = {k: list(v) for k, v in pathway_to_genes.items()}
    adata.uns["source"] = "msep.datasets.load_example"
    adata.uns["seed"] = seed
    return adata


def list_examples() -> Dict[str, str]:
    """Available built-in example datasets.

    Returns a mapping from dataset name to a one-line description. Only
    one dataset is bundled in msep 1.x; more may be added in future
    releases.
    """
    return {
        "chordoma_demo": (
            "Synthetic 500-cell chordoma-like dataset reproducing the "
            "multi-scale paradox (high CSC entropy + low EMT CV + high "
            "immune-evasion CV)."
        ),
    }
