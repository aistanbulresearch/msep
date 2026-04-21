"""
Result container for multi-scale entropy profiling.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


@dataclass
class MSEPResult:
    """Container returned by :func:`msep.profile`.

    Attributes
    ----------
    per_cell_entropy : DataFrame
        Per-cell entropy values with columns for global and pathway-
        specific entropy, plus normalized versions.
    pathway_cv : DataFrame
        Across-cells CV per (cell_type, pathway).
    gene_cv : dict[str, DataFrame] or None
        Per-gene CV tables, keyed by pathway name.
    bootstrap : dict[str, dict] or None
        Bootstrap CI results, keyed by ``"cell_type|pathway"``.
    perturbation : DataFrame or None
        Pseudo-perturbation results.
    xbp1 : DataFrame or None
        XBP1 consolidation results.
    pathways_used : dict[str, list[str]]
        The pathway gene sets that were actually used.
    metadata : dict
        Analysis parameters and summary statistics.
    """

    per_cell_entropy: pd.DataFrame
    pathway_cv: pd.DataFrame
    gene_cv: Optional[Dict[str, pd.DataFrame]] = None
    bootstrap: Optional[Dict[str, dict]] = None
    perturbation: Optional[pd.DataFrame] = None
    xbp1: Optional[pd.DataFrame] = None
    pathways_used: Dict[str, List[str]] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    # -------------------------------------------------------------------
    # Convenience accessors
    # -------------------------------------------------------------------

    @property
    def coordination_scores(self) -> pd.DataFrame:
        """Per ``(cell_type, pathway)`` coordination scores with classification.

        Wraps :func:`msep.scoring.coordination_score_table` with the
        current result's entropy and CV tables, then adds a
        ``classification`` column via :func:`msep.scoring.classify_paradox`.

        Returns an empty DataFrame if the result is missing the underlying
        tables (e.g. after a partial pipeline run).
        """
        from .scoring import coordination_score_table, classify_paradox

        if self.per_cell_entropy.empty or self.pathway_cv.empty:
            return pd.DataFrame(
                columns=["cell_type", "pathway", "coordination_score",
                         "composite_score", "classification"]
            )

        table = coordination_score_table(self.per_cell_entropy,
                                         self.pathway_cv)
        if table.empty:
            return table
        table = table.copy()
        table["classification"] = table["coordination_score"].map(classify_paradox)
        return table

    @property
    def paradox_summary(self) -> pd.DataFrame:
        """One-row-per-cell-type summary: per-cell entropy vs. across-cells CV.

        Useful for identifying the 'individually diverse, collectively
        disciplined' phenotype.
        """
        if self.per_cell_entropy.empty or self.pathway_cv.empty:
            return pd.DataFrame()

        entropy_medians = (
            self.per_cell_entropy
            .groupby("cell_type")["entropy_global"]
            .median()
            .rename("median_entropy")
        )

        cv_pivot = self.pathway_cv.pivot(
            index="cell_type", columns="pathway", values="cv"
        )
        cv_pivot.columns = [f"cv_{c}" for c in cv_pivot.columns]

        return entropy_medians.to_frame().join(cv_pivot)

    def consolidation_score(self, baseline: str = "XBP1-zero") -> Optional[pd.DataFrame]:
        """Number of pathways showing consolidation (lower CV) under stress.

        Returns
        -------
        DataFrame or None
            Columns: ``pathway``, ``cv_baseline``, ``cv_high``,
            ``delta``, ``consolidated`` (bool).
        """
        if self.xbp1 is None or self.xbp1.empty:
            return None

        base = self.xbp1[self.xbp1["xbp1_group"] == baseline].set_index("pathway")["cv"]
        high = self.xbp1[self.xbp1["xbp1_group"] == "XBP1-high"].set_index("pathway")["cv"]

        merged = pd.DataFrame({
            "cv_baseline": base,
            "cv_high": high,
        }).dropna()
        merged["delta"] = merged["cv_high"] - merged["cv_baseline"]
        merged["consolidated"] = merged["delta"] < 0
        return merged.reset_index()

    # -------------------------------------------------------------------
    # Repr
    # -------------------------------------------------------------------

    def __repr__(self) -> str:
        n_cells = len(self.per_cell_entropy)
        n_ct = self.per_cell_entropy["cell_type"].nunique() if n_cells else 0
        n_pw = self.pathway_cv["pathway"].nunique() if not self.pathway_cv.empty else 0
        parts = [
            f"MSEPResult(n_cells={n_cells:,}, cell_types={n_ct}, pathways={n_pw}",
        ]
        if self.perturbation is not None:
            parts.append(f"  perturbation={len(self.perturbation)} tests")
        if self.xbp1 is not None:
            parts.append(f"  xbp1_groups={self.xbp1['xbp1_group'].nunique()}")
        parts.append(")")
        return "\n".join(parts)
