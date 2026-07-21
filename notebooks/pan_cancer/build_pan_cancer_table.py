"""Derive the committed pan-cancer CV table from the analysis output.

The heavy pipeline (CellxGene Census queries, GEO downloads, malignant-cell
selection, ``msep`` profiling) is documented in ``examples/chordoma_msep.ipynb``
and writes its raw output to ``processed/``. That directory is large and is not
tracked in git.

This script distils the one table Figure 3 needs into a small, tracked,
provenance-carrying CSV at ``data/pan_cancer_cv.csv`` so that
``notebooks/figures/figure3_pan_cancer.ipynb`` can render the figure in seconds
without re-running the pipeline.

Usage::

    python notebooks/pan_cancer/build_pan_cancer_table.py

Re-run this after extending the pan-cancer panel (see WP-3 in
``plan/GB_SUBMISSION_WORKPLAN.md``); the figure and its tests pick up the new
rows automatically.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_CSV = REPO_ROOT / "processed" / "pan_cancer_census_cv_expanded.csv"
OUTPUT_CSV = REPO_ROOT / "data" / "pan_cancer_cv.csv"

# Provenance for each row of the pan-cancer panel. The Census block was queried
# from the CellxGene Census stable release; the remaining three came from GEO.
# Keyed by the ``cancer_type`` value used in the pipeline output.
PROVENANCE: dict[str, dict[str, object]] = {
    "Chordoma": {
        "display_name": "Chordoma",
        "source": "This study",
        "accession": "Zenodo (Arrieta et al. 2025)",
        "platform": "10x",
        "cell_selection": "CSC_TBXT+",
    },
    "GBM": {
        "display_name": "Glioblastoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Lung_Adeno": {
        "display_name": "Lung adenocarcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Breast_Ductal": {
        "display_name": "Invasive ductal breast carcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Breast_TNBC": {
        "display_name": "Triple-negative breast carcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "RCC_Nonpap": {
        "display_name": "Nonpapillary renal cell carcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Lung_SCC": {
        "display_name": "Squamous cell lung carcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Melanoma_Uveal": {
        "display_name": "Uveal melanoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "BCC": {
        "display_name": "Basal cell carcinoma",
        "source": "CellxGene Census",
        "accession": "census_2025-11-08",
        "platform": "10x",
        "cell_selection": "malignant (author-labeled)",
    },
    "Melanoma": {
        "display_name": "Melanoma",
        "source": "GEO",
        "accession": "GSE115978",
        "platform": "Smart-seq2",
        "cell_selection": "Mal (annotated)",
    },
    "Breast": {
        "display_name": "Breast cancer",
        "source": "GEO",
        "accession": "GSE176078",
        "platform": "10x",
        "cell_selection": "Cancer Epithelial",
    },
    "Osteosarcoma*": {
        "display_name": "Osteosarcoma",
        "source": "GEO",
        "accession": "GSE162454",
        "platform": "10x",
        "cell_selection": "All cells (no annotation)",
    },
}

# ``emt_mechanical`` in the pipeline output is the pathway the manuscript and the
# ``msep`` built-in gene sets both call ``emt``.
PATHWAY_RENAMES = {
    "ferroptosis_CV": "ferroptosis_cv",
    "ferroptosis_n_genes": "ferroptosis_n_genes",
    "immune_evasion_CV": "immune_evasion_cv",
    "immune_evasion_n_genes": "immune_evasion_n_genes",
    "emt_mechanical_CV": "emt_cv",
    "emt_mechanical_n_genes": "emt_n_genes",
    "housekeeping_CV": "housekeeping_cv",
    "housekeeping_n_genes": "housekeeping_n_genes",
}

COLUMN_ORDER = [
    "cancer_type",
    "display_name",
    "source",
    "accession",
    "platform",
    "n_cells",
    "cell_selection",
    "in_primary_comparison",
    "emt_cv",
    "emt_n_genes",
    "ferroptosis_cv",
    "ferroptosis_n_genes",
    "immune_evasion_cv",
    "immune_evasion_n_genes",
    "housekeeping_cv",
    "housekeeping_n_genes",
]


def build_table(source_csv: Path = SOURCE_CSV) -> pd.DataFrame:
    """Read the pipeline output and return the tidy, provenance-carrying table.

    Args:
        source_csv: Path to ``pan_cancer_census_cv_expanded.csv``.

    Returns:
        The pan-cancer table ready to be written to ``data/pan_cancer_cv.csv``.

    Raises:
        FileNotFoundError: If the pipeline output is missing.
        ValueError: If a cancer type has no registered provenance.
    """
    if not source_csv.exists():
        raise FileNotFoundError(
            f"Pipeline output not found: {source_csv}. "
            "Run the pan-cancer section of examples/chordoma_msep.ipynb first."
        )

    df = pd.read_csv(source_csv).rename(columns=PATHWAY_RENAMES)

    unknown = set(df["cancer_type"]) - set(PROVENANCE)
    if unknown:
        raise ValueError(
            f"No provenance registered for {sorted(unknown)}. "
            "Add an entry to PROVENANCE before rebuilding the table."
        )

    meta = pd.DataFrame.from_dict(PROVENANCE, orient="index").rename_axis("cancer_type")
    # Provenance is authoritative: the pipeline output leaves source/platform
    # blank on the three GEO rows that were merged in from an earlier run.
    df = df.drop(columns=["source", "platform", "cell_selection", "disease"], errors="ignore")
    df = df.merge(meta, on="cancer_type", how="left")

    # The manuscript restricts the primary pan-cancer comparison to 10x data;
    # Smart-seq2 melanoma is reported separately because of known platform
    # effects on expression variance (Methods 2.1.3, Results 3.4).
    df["in_primary_comparison"] = df["platform"].eq("10x")

    df = df[COLUMN_ORDER].sort_values("emt_cv").reset_index(drop=True)
    return df


def main() -> None:
    """Build the table and write it to ``data/pan_cancer_cv.csv``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    df = build_table()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    logger.info("Wrote %d cancer types to %s", len(df), OUTPUT_CSV)
    logger.info(
        "  %d in primary (10x) comparison, %d reported separately",
        int(df["in_primary_comparison"].sum()),
        int((~df["in_primary_comparison"]).sum()),
    )


if __name__ == "__main__":
    main()
