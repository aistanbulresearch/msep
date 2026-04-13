"""
Per-cell Shannon entropy computation on raw count matrices.

The probability of observing a read from gene *i* in a given cell is

    p_i = count_i / total_counts

and Shannon entropy is

    H = −Σ p_i log₂(p_i)   (bits)

Genes with zero counts are excluded.  The calculation is performed on
**raw (unnormalized)** counts to avoid normalization-induced artifacts.
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix

__all__ = ["per_cell_entropy"]


def per_cell_entropy(
    count_matrix: Union[np.ndarray, spmatrix],
    gene_indices: Optional[np.ndarray] = None,
    exclude_indices: Optional[np.ndarray] = None,
    base: int = 2,
    min_total_counts: int = 10,
    chunk_size: int = 5000,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-cell Shannon entropy from a raw count matrix.

    Parameters
    ----------
    count_matrix
        Cells × genes matrix of **raw** (unnormalized) counts.
        Sparse (CSR/CSC) or dense.
    gene_indices
        If given, restrict computation to these column indices.
    exclude_indices
        If given (and *gene_indices* is ``None``), exclude these columns
        (e.g. ribosomal genes).
    base
        Logarithm base.  ``2`` → bits (default), ``e`` → nats.
    min_total_counts
        Cells with total counts below this threshold receive ``NaN``.
    chunk_size
        Number of cells processed per chunk (memory control).

    Returns
    -------
    entropy : ndarray, shape (n_cells,)
        Shannon entropy per cell.
    n_expressed : ndarray, shape (n_cells,)
        Number of genes with count > 0 per cell.
    """
    n_cells = count_matrix.shape[0]

    # --- gene subsetting ------------------------------------------------
    mat = _subset_genes(count_matrix, gene_indices, exclude_indices)

    # --- allocate output -------------------------------------------------
    entropy = np.full(n_cells, np.nan, dtype=np.float64)
    n_expressed = np.zeros(n_cells, dtype=np.int32)
    log_fn = np.log2 if base == 2 else np.log

    # --- chunked computation --------------------------------------------
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = mat[start:end]
        if issparse(chunk):
            chunk = chunk.toarray()
        chunk = chunk.astype(np.float64, copy=False)

        totals = chunk.sum(axis=1)

        for i in range(chunk.shape[0]):
            idx = start + i
            total = totals[i]
            if total < min_total_counts:
                continue

            row = chunk[i]
            mask = row > 0
            n_expr = int(mask.sum())
            n_expressed[idx] = n_expr

            if n_expr <= 1:
                entropy[idx] = 0.0
                continue

            p = row[mask] / total
            entropy[idx] = -np.dot(p, log_fn(p))

    return entropy, n_expressed


def normalized_entropy(
    entropy: np.ndarray,
    n_expressed: np.ndarray,
    base: int = 2,
) -> np.ndarray:
    """Normalize entropy to [0, 1] by dividing by log(n_expressed).

    Parameters
    ----------
    entropy
        Raw entropy values (from :func:`per_cell_entropy`).
    n_expressed
        Number of expressed genes per cell.
    base
        Must match the *base* used in :func:`per_cell_entropy`.

    Returns
    -------
    ndarray
        Normalized entropy.  ``NaN`` where *n_expressed* ≤ 1.
    """
    log_fn = np.log2 if base == 2 else np.log
    max_h = np.where(n_expressed > 1, log_fn(n_expressed.astype(float)), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.isnan(max_h), np.nan, entropy / max_h)


# -----------------------------------------------------------------------
# internal helpers
# -----------------------------------------------------------------------

def _subset_genes(
    mat: Union[np.ndarray, spmatrix],
    gene_indices: Optional[np.ndarray],
    exclude_indices: Optional[np.ndarray],
) -> Union[np.ndarray, spmatrix]:
    if gene_indices is not None:
        return mat[:, gene_indices]
    if exclude_indices is not None:
        keep = np.setdiff1d(np.arange(mat.shape[1]), exclude_indices)
        return mat[:, keep]
    return mat
