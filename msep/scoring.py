"""
Coordination Score: a single metric summarising the multi-scale paradox.

The MSEP Coordination Score quantifies the degree to which a cell
population is "individually diverse, collectively disciplined" for a
given pathway or set of pathways.

Definition
----------
For a cell type *c* and pathway *p*:

    CS(c, p) = entropy_rank(c) × (1 − cv_rank(c, p))

where:
- ``entropy_rank`` is the fractional rank of *c*'s median per-cell
  entropy among all cell types (1.0 = highest entropy)
- ``cv_rank`` is the fractional rank of *c*'s pathway CV among all
  cell types (1.0 = highest CV, i.e. least coordinated)

The score ranges from 0 to 1:
- **CS → 1**: high per-cell entropy + low across-cells CV
  = individually diverse, collectively disciplined (paradox phenotype)
- **CS → 0**: low entropy and/or high CV = no paradox

The composite score across all pathways is the mean of per-pathway
scores, optionally weighted.

This metric enables direct, quantitative comparison of the paradox
phenotype across cell types and across cancer types.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

__all__ = ["coordination_score", "coordination_score_table"]


def coordination_score(
    entropy_medians: pd.Series,
    pathway_cvs: pd.Series,
) -> pd.DataFrame:
    """Compute coordination score for each cell type on one pathway.

    Parameters
    ----------
    entropy_medians
        Median per-cell entropy per cell type.  Index = cell type names.
    pathway_cvs
        Pathway CV per cell type.  Index = cell type names.

    Returns
    -------
    DataFrame
        Columns: ``cell_type``, ``entropy_rank``, ``cv_rank``,
        ``coordination_score``.
    """
    common = entropy_medians.index.intersection(pathway_cvs.index)
    if len(common) < 2:
        return pd.DataFrame(
            columns=["cell_type", "entropy_rank", "cv_rank",
                     "coordination_score"]
        )

    e = entropy_medians.loc[common]
    c = pathway_cvs.loc[common]
    n = len(common)

    # Fractional ranks (higher entropy → higher rank)
    e_rank = e.rank(method="average") / n
    # Higher CV → higher rank (less coordinated)
    c_rank = c.rank(method="average") / n

    # CS = entropy_rank × (1 - cv_rank)
    # High entropy + low CV → high score
    cs = e_rank * (1.0 - c_rank)

    return pd.DataFrame({
        "cell_type": common,
        "median_entropy": e.values,
        "pathway_cv": c.values,
        "entropy_rank": np.round(e_rank.values, 4),
        "cv_rank": np.round(c_rank.values, 4),
        "coordination_score": np.round(cs.values, 4),
    }).sort_values("coordination_score", ascending=False).reset_index(drop=True)


def coordination_score_table(
    per_cell_entropy: pd.DataFrame,
    pathway_cv: pd.DataFrame,
    pathways: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Compute coordination score for all (cell_type, pathway) combinations.

    Parameters
    ----------
    per_cell_entropy
        DataFrame with columns ``cell_type`` and ``entropy_global``
        (as returned by :func:`msep.profile`).
    pathway_cv
        DataFrame with columns ``cell_type``, ``pathway``, ``cv``
        (as returned by :func:`msep.profile`).
    pathways
        Subset of pathways to score.  If ``None``, scores all.

    Returns
    -------
    DataFrame
        Columns: ``cell_type``, ``pathway``, ``entropy_rank``,
        ``cv_rank``, ``coordination_score``, ``composite_score``.
        The ``composite_score`` is the mean coordination score across
        all pathways for each cell type.
    """
    entropy_medians = (
        per_cell_entropy
        .groupby("cell_type")["entropy_global"]
        .median()
    )

    available_pathways = pathway_cv["pathway"].unique()
    if pathways is not None:
        available_pathways = [p for p in pathways if p in available_pathways]

    rows = []
    for pw in available_pathways:
        pw_data = pathway_cv[pathway_cv["pathway"] == pw]
        pw_cvs = pw_data.set_index("cell_type")["cv"].dropna()

        if len(pw_cvs) < 2:
            continue

        cs_df = coordination_score(entropy_medians, pw_cvs)
        cs_df["pathway"] = pw
        rows.append(cs_df)

    if not rows:
        return pd.DataFrame(
            columns=["cell_type", "pathway", "median_entropy",
                     "pathway_cv", "entropy_rank", "cv_rank",
                     "coordination_score", "composite_score"]
        )

    result = pd.concat(rows, ignore_index=True)

    # Composite score: mean across pathways per cell type
    composite = (
        result
        .groupby("cell_type")["coordination_score"]
        .mean()
        .rename("composite_score")
    )
    result = result.merge(composite, on="cell_type")

    # Sort by composite then per-pathway score
    result = result.sort_values(
        ["composite_score", "coordination_score"],
        ascending=[False, False],
    ).reset_index(drop=True)

    return result


def classify_paradox(
    score: float,
    high_threshold: float = 0.7,
    moderate_threshold: float = 0.4,
) -> str:
    """Classify a coordination score into a paradox category.

    Parameters
    ----------
    score
        Coordination score (0–1).
    high_threshold
        Above this → "strong paradox".
    moderate_threshold
        Above this → "moderate paradox".

    Returns
    -------
    str
        One of: ``"strong paradox"``, ``"moderate paradox"``,
        ``"weak/absent"``.
    """
    if np.isnan(score):
        return "insufficient data"
    if score >= high_threshold:
        return "strong paradox"
    elif score >= moderate_threshold:
        return "moderate paradox"
    else:
        return "weak/absent"
