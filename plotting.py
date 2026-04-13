"""
Publication-ready plotting for MSEP results.

All functions accept an :class:`~msep.result.MSEPResult` and return a
``matplotlib.figure.Figure``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import matplotlib.figure

__all__ = [
    "plot_entropy_violin",
    "plot_pathway_cv_heatmap",
    "plot_paradox",
    "plot_pan_cancer",
]

# -----------------------------------------------------------------------
# style defaults
# -----------------------------------------------------------------------

PATHWAY_COLORS = {
    "ferroptosis": "#E74C3C",
    "immune_evasion": "#3498DB",
    "emt": "#2ECC71",
    "housekeeping": "#95A5A6",
}

_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.dpi": 150,
}


def _apply_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update(_RC)


# -----------------------------------------------------------------------
# 1. Per-cell entropy violin
# -----------------------------------------------------------------------

def plot_entropy_violin(
    result,
    key: str = "entropy_global",
    order: Optional[List[str]] = None,
    figsize: tuple = (12, 5),
):
    """Violin plot of per-cell entropy across cell types.

    Parameters
    ----------
    result : MSEPResult
    key : str
        Column in ``result.per_cell_entropy`` to plot.
    order : list[str] or None
        Cell-type order (default: descending median).
    figsize : tuple

    Returns
    -------
    Figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    _apply_style()

    df = result.per_cell_entropy[["cell_type", key]].dropna()
    if order is None:
        order = (
            df.groupby("cell_type")[key]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )

    fig, ax = plt.subplots(figsize=figsize)
    sns.violinplot(
        data=df, x="cell_type", y=key, order=order,
        inner="quartile", cut=0, ax=ax, palette="Set2",
    )
    ax.set_title(f"Per-cell entropy by cell type ({key})")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=45)

    # annotate medians
    for i, ct in enumerate(order):
        med = df.loc[df["cell_type"] == ct, key].median()
        ax.text(i, med + 0.1, f"{med:.2f}", ha="center", fontsize=8,
                fontweight="bold")

    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------
# 2. CV heatmap
# -----------------------------------------------------------------------

def plot_pathway_cv_heatmap(
    result,
    cell_types: Optional[List[str]] = None,
    figsize: tuple = (10, 5),
    cmap: str = "YlOrRd_r",
    vmax: float = 12,
):
    """Heatmap of across-cells CV: cell type × pathway.

    Parameters
    ----------
    result : MSEPResult
    cell_types : list[str] or None
    figsize, cmap, vmax : display settings

    Returns
    -------
    Figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    _apply_style()

    df = result.pathway_cv.copy()
    if cell_types is not None:
        df = df[df["cell_type"].isin(cell_types)]

    pivot = df.pivot(index="cell_type", columns="pathway", values="cv")
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot, annot=True, fmt=".2f", cmap=cmap,
        linewidths=0.5, vmin=0, vmax=vmax, ax=ax,
        cbar_kws={"label": "CV (lower = more coordinated)"},
    )
    ax.set_title("Across-cells CV: cell type × pathway")
    ax.set_ylabel("")
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------
# 3. Paradox summary (entropy vs CV scatter)
# -----------------------------------------------------------------------

def plot_paradox(
    result,
    cv_pathway: str = "emt",
    figsize: tuple = (8, 6),
):
    """Scatter plot of per-cell entropy (median) vs across-cells CV.

    Each point is a cell type.  The 'individually diverse, collectively
    disciplined' phenotype occupies the top-left quadrant (high entropy,
    low CV).

    Parameters
    ----------
    result : MSEPResult
    cv_pathway : str
        Pathway to plot on the y-axis.
    figsize : tuple

    Returns
    -------
    Figure
    """
    import matplotlib.pyplot as plt
    _apply_style()

    summary = result.paradox_summary
    cv_col = f"cv_{cv_pathway}"
    if cv_col not in summary.columns:
        raise ValueError(
            f"No CV column {cv_col!r}. Available: {list(summary.columns)}"
        )

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        summary["median_entropy"], summary[cv_col],
        s=120, edgecolors="black", linewidths=0.8,
        c=range(len(summary)), cmap="tab10", zorder=3,
    )
    for i, ct in enumerate(summary.index):
        ax.annotate(
            ct,
            (summary.loc[ct, "median_entropy"], summary.loc[ct, cv_col]),
            textcoords="offset points", xytext=(6, 6), fontsize=9,
        )

    ax.set_xlabel("Median per-cell entropy (bits)")
    ax.set_ylabel(f"Across-cells CV ({cv_pathway})")
    ax.set_title(
        "Multi-scale paradox: per-cell diversity vs population coordination"
    )
    ax.axhline(summary[cv_col].median(), ls="--", color="gray", alpha=0.4)
    ax.axvline(summary["median_entropy"].median(), ls="--", color="gray", alpha=0.4)

    # Quadrant label
    ax.text(
        0.02, 0.02,
        "← Individually diverse,\n    collectively disciplined",
        transform=ax.transAxes, fontsize=9, fontstyle="italic",
        color="#8E44AD",
    )
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------
# 4. Pan-cancer bar chart
# -----------------------------------------------------------------------

def plot_pan_cancer(
    cv_data: pd.DataFrame,
    pathway: str = "emt",
    highlight: Optional[str] = None,
    figsize: tuple = (10, 5),
):
    """Horizontal bar chart ranking cancer types by pathway CV.

    Parameters
    ----------
    cv_data : DataFrame
        Must have columns ``cancer_type`` and ``cv``.
    pathway : str
        Label for the title.
    highlight : str or None
        Cancer type name to highlight in a different colour.
    figsize : tuple

    Returns
    -------
    Figure
    """
    import matplotlib.pyplot as plt
    _apply_style()

    df = cv_data.sort_values("cv", ascending=True).reset_index(drop=True)
    colors = [
        "#E74C3C" if ct == highlight else "#3498DB"
        for ct in df["cancer_type"]
    ]

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(df["cancer_type"], df["cv"], color=colors, edgecolor="white")
    ax.set_xlabel(f"Across-cells CV ({pathway})")
    ax.set_title(f"{pathway.upper()} population coordination ranking")

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["cv"] + 0.05, i, f"{row['cv']:.2f}", va="center", fontsize=9)

    fig.tight_layout()
    return fig
