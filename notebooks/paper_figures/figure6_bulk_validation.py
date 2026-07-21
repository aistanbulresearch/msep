"""Regenerate manuscript Figure 6 (bulk RNA-seq validation) from analysis output.

Every number printed on these figures is computed from the CSVs in ``processed/``
rather than typed in.

The motivating defect was the old panel D. The version carried in the manuscript
reported ``Spearman rho = 0.000, p = 1.0000`` over an 8-gene subset. That statistic
matched neither the source data (23 genes, rho = -0.110, p = 0.618) nor the eight
points the panel itself plotted (rho = -0.500, p = 0.207), and several plotted
coordinates were wrong -- SLC7A11 appeared at a scRNA-seq CV of ~2.1 against a true
value of 13.73.

Recomputed honestly over all 23 genes the correlation is weak and not significant,
so gene-level CV does not concord across platforms. The claim section 3.7 actually
rests on is the *pathway*-level hierarchy, which panel A establishes and which is
unaffected. Figure 6 therefore keeps three panels, and the gene-level scatter moves
to the supplement, where a null result can be reported plainly.

Usage::

    python notebooks/paper_figures/figure6_bulk_validation.py

Writes:
    figures/Publication/Figure_6_BulkValidation.{pdf,png}          (panels A-C)
    figures/Supplementary/Figure_S15_GeneLevel_CV_Concordance.{pdf,png}

``Figure_6D_GeneLevel_CV_Concordance_REAL.{pdf,png}`` in ``figures/Supplementary/``
is superseded by the S15 output and can be removed.
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
PUBLICATION_DIR = REPO_ROOT / "figures" / "Publication"
SUPPLEMENTARY_DIR = REPO_ROOT / "figures" / "Supplementary"

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


def supplementary_gene_level(ax: plt.Axes) -> tuple[float, float]:
    """Gene-level CV, scRNA-seq vs bulk, over all 23 genes.

    Returns:
        The Spearman correlation and p-value, so the caller can log them.
    """
    genes = _read("gene_cv_concordance_scRNA_vs_bulk.csv")
    rho, pval = spearmanr(genes["scRNA_CV"], genes["Bulk_CV"])

    for category, sub in genes.groupby("category"):
        ax.scatter(sub["scRNA_CV"], sub["Bulk_CV"], s=55,
                   color=CATEGORY_COLORS.get(category, '#7F8C8D'),
                   edgecolor='black', linewidth=0.6, zorder=3,
                   label=category.replace('_', ' '))

    # Three pairs sit close enough that their default labels overlap; nudge one of
    # each pair rather than pulling in a label-layout dependency.
    label_offsets = {
        "CD47": (5, -9),
        "SMAD2": (5, -9),
        "SMAD3": (5, 4),
        "TGFBR2": (-38, -6),
    }
    for _, row in genes.iterrows():
        ax.annotate(row["gene"], (row["scRNA_CV"], row["Bulk_CV"]),
                    textcoords='offset points',
                    xytext=label_offsets.get(row["gene"], (5, 3)), fontsize=7)

    # SLC7A11 and CD274 sit an order of magnitude out on the scRNA-seq axis; a log
    # scale keeps the remaining 21 genes legible instead of collapsing them.
    ax.set_xscale('log')
    ax.set_xlabel('scRNA-seq gene CV (CSC_TBXT+, across cells), log scale')
    ax.set_ylabel('Bulk RNA-seq gene CV (across tumors)')

    verdict = 'not significant' if pval >= 0.05 else 'significant'
    ax.set_title(
        f'Gene-level CV concordance, scRNA-seq vs bulk RNA-seq\n'
        f'Spearman $\\rho$ = {rho:.3f}, p = {pval:.3f} ({verdict}), n = {len(genes)}'
    )
    ax.legend(fontsize=8, loc='best')
    ax.spines[['top', 'right']].set_visible(False)
    return float(rho), float(pval)


def build_figure6() -> None:
    """Panels A-C: the pathway-level validation Figure 6 now carries."""
    # A and B share the pathway axis and sit side by side; C spans the width below
    # because its six patient labels need the room.
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[1, :])]

    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])

    for ax, letter in zip(axes, "ABC"):
        ax.text(-0.06, 1.06, letter, transform=ax.transAxes,
                fontsize=17, fontweight='bold', va='top')

    fig.tight_layout()
    PUBLICATION_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PUBLICATION_DIR / "Figure_6_BulkValidation.pdf", bbox_inches='tight')
    fig.savefig(PUBLICATION_DIR / "Figure_6_BulkValidation.png", dpi=300,
                bbox_inches='tight')
    plt.close(fig)
    logger.info("Wrote Figure_6_BulkValidation.{pdf,png} to %s", PUBLICATION_DIR)


def build_supplementary() -> tuple[float, float]:
    """Figure S15: the gene-level scatter, reported as the null result it is."""
    fig, ax = plt.subplots(figsize=(9, 7))
    rho, pval = supplementary_gene_level(ax)

    fig.tight_layout()
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)
    stem = SUPPLEMENTARY_DIR / "Figure_S15_GeneLevel_CV_Concordance"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches='tight')
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("Wrote Figure_S15_GeneLevel_CV_Concordance.{pdf,png} to %s",
                SUPPLEMENTARY_DIR)
    return rho, pval


def main() -> None:
    """Build Figure 6 (A-C) and Supplementary Figure S15."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    build_figure6()
    rho, pval = build_supplementary()
    logger.info("Gene-level Spearman rho = %.3f, p = %.3f", rho, pval)


if __name__ == "__main__":
    main()
