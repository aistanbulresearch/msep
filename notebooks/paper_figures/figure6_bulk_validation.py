"""Regenerate manuscript Figure 6 (bulk RNA-seq validation) from analysis output.

This rebuilds all four panels from the CSVs in ``processed/`` so that every number
printed on the figure is derived, not typed.

The motivating defect was panel D. The version carried in the manuscript reported
``Spearman rho = 0.000, p = 1.0000`` over an 8-gene subset. That statistic matched
neither the source data (23 genes, rho = -0.110, p = 0.618) nor the eight points the
panel itself plotted (rho = -0.500, p = 0.207), and several plotted coordinates were
wrong -- SLC7A11 appeared at a scRNA-seq CV of ~2.1 against a true value of 13.73.

Panel D is now computed from all 23 genes in
``processed/gene_cv_concordance_scRNA_vs_bulk.csv`` and reports the correlation it
actually finds, including when that correlation is not significant.

Usage::

    python notebooks/paper_figures/figure6_bulk_validation.py

Writes ``figures/Publication/Figure_6_BulkValidation.{pdf,png}``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
PROCESSED = REPO_ROOT / "processed"
OUTPUT_DIR = REPO_ROOT / "figures" / "Publication"

# Pathway display order, shared across panels A and B.
PATHWAY_ORDER = ["emt_mechanical", "housekeeping", "ferroptosis", "immune_evasion"]
PATHWAY_LABELS = ["EMT", "Housekeeping", "Ferroptosis", "Immune\nEvasion"]

# Across-cells CV for chordoma CSC_TBXT+, the scRNA-seq side of panel A.
SCRNA_CV_COLUMN = "cancer_type"

CATEGORY_COLORS = {
    "EMT_shield": "#E74C3C",
    "Ferroptosis_shield": "#F39C12",
    "Immune_shield": "#3498DB",
    "TGFb_axis": "#9B59B6",
}


def _read(name: str) -> pd.DataFrame:
    """Read a CSV from ``processed/``, with a pointed error if it is absent."""
    path = PROCESSED / name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. This script regenerates the manuscript figure from "
            "the full analysis output, which is not tracked in git."
        )
    return pd.read_csv(path)


def panel_a(ax: plt.Axes) -> None:
    """Cross-platform pathway CV concordance: scRNA-seq vs bulk."""
    bulk = _read("bulk_validation_cv.csv").set_index("group")
    pan = _read("pan_cancer_census_cv_expanded.csv").set_index("cancer_type")

    scrna = [pan.loc["Chordoma", f"{p}_CV"] for p in PATHWAY_ORDER]
    bulk_cv = [bulk.loc["Chordoma", f"{p}_CV"] for p in PATHWAY_ORDER]

    x = range(len(PATHWAY_ORDER))
    w = 0.38
    ax.bar([i - w/2 for i in x], scrna, w, label='scRNA-seq (across-cells)',
           color='#E74C3C', edgecolor='black', linewidth=0.6)
    ax.bar([i + w/2 for i in x], bulk_cv, w, label='Bulk RNA-seq (across-tumors)',
           color='#3498DB', edgecolor='black', linewidth=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels(PATHWAY_LABELS)
    ax.set_ylabel('CV')
    ax.set_title('Cross-platform pathway CV concordance')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def panel_b(ax: plt.Axes) -> None:
    """Chordoma vs notochord bulk CV -- the developmental-origin comparison."""
    bulk = _read("bulk_validation_cv.csv").set_index("group")

    chordoma = [bulk.loc["Chordoma", f"{p}_CV"] for p in PATHWAY_ORDER]
    notochord = [bulk.loc["Notochord", f"{p}_CV"] for p in PATHWAY_ORDER]

    x = range(len(PATHWAY_ORDER))
    w = 0.38
    ax.bar([i - w/2 for i in x], chordoma, w, label='Chordoma (n=6)',
           color='#E74C3C', edgecolor='black', linewidth=0.6)
    ax.bar([i + w/2 for i in x], notochord, w, label='Notochord (n=2)',
           color='#2ECC71', edgecolor='black', linewidth=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels(PATHWAY_LABELS)
    ax.set_ylabel('Bulk CV (across tumors/samples)')
    ax.set_title('Chordoma vs Notochord (developmental origin)')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def panel_c(ax: plt.Axes) -> None:
    """Per-patient pathway CV stability, starred where EMT is most homogeneous."""
    per_patient = _read("per_patient_pathway_cv.csv")

    resistance = ["emt_mechanical", "ferroptosis", "immune_evasion"]
    labels = ["EMT", "Ferroptosis", "Immune"]
    colors = ["#E74C3C", "#F39C12", "#3498DB"]

    x = range(len(per_patient))
    w = 0.26
    for offset, (pw, label, color) in enumerate(zip(resistance, labels, colors)):
        ax.bar([i + (offset - 1) * w for i in x], per_patient[f"{pw}_CV"], w,
               label=label, color=color, edgecolor='black', linewidth=0.5)

    # Star the patients where EMT is the most homogeneous resistance pathway.
    for i, row in per_patient.iterrows():
        if row["most_homogeneous"] == "emt_mechanical":
            top = max(row[f"{pw}_CV"] for pw in resistance)
            ax.text(i, top + 0.25, '*', ha='center', fontsize=15, fontweight='bold')

    ax.set_xticks(list(x))
    ax.set_xticklabels(per_patient["patient_id"], rotation=45, ha='right')
    ax.set_ylabel('Across-cells CV')
    ax.set_title('Per-patient pathway CV stability')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)


def panel_d(ax: plt.Axes) -> float:
    """Gene-level CV, scRNA-seq vs bulk, over all 23 genes.

    Returns:
        The Spearman correlation actually computed, so the caller can log it.
    """
    genes = _read("gene_cv_concordance_scRNA_vs_bulk.csv")
    rho, pval = spearmanr(genes["scRNA_CV"], genes["Bulk_CV"])

    for category, sub in genes.groupby("category"):
        ax.scatter(sub["scRNA_CV"], sub["Bulk_CV"], s=55,
                   color=CATEGORY_COLORS.get(category, '#7F8C8D'),
                   edgecolor='black', linewidth=0.6, zorder=3,
                   label=category.replace('_', ' '))

    for _, row in genes.iterrows():
        ax.annotate(row["gene"], (row["scRNA_CV"], row["Bulk_CV"]),
                    textcoords='offset points', xytext=(5, 3), fontsize=6.5)

    # SLC7A11 and CD274 sit an order of magnitude out on the scRNA-seq axis; a log
    # scale keeps the remaining 21 genes legible instead of collapsing them.
    ax.set_xscale('log')
    ax.set_xlabel('scRNA-seq gene CV (CSC_TBXT+, across cells), log scale')
    ax.set_ylabel('Bulk RNA-seq gene CV (across tumors)')

    significance = 'n.s.' if pval >= 0.05 else f'p = {pval:.3f}'
    ax.set_title(
        f'Gene-level CV concordance\n'
        f'(Spearman $\\rho$ = {rho:.3f}, p = {pval:.3f}, {significance}, n = {len(genes)})'
    )
    ax.legend(fontsize=7, loc='best')
    ax.spines[['top', 'right']].set_visible(False)
    return float(rho)


def main() -> None:
    """Build all four panels and write the figure."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[1, 0])
    rho = panel_d(axes[1, 1])

    for ax, letter in zip(axes.flat, "ABCD"):
        ax.text(-0.09, 1.06, letter, transform=ax.transAxes,
                fontsize=17, fontweight='bold', va='top')

    fig.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / "Figure_6_BulkValidation.pdf", bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "Figure_6_BulkValidation.png", dpi=300, bbox_inches='tight')

    logger.info("Wrote Figure_6_BulkValidation.{pdf,png} to %s", OUTPUT_DIR)
    logger.info("Panel D gene-level Spearman rho = %.3f", rho)


if __name__ == "__main__":
    main()
