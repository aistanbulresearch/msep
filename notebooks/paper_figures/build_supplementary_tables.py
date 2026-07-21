"""Build Supplementary Tables S1-S10 for the manuscript.

The manuscript declares ten supplementary tables but none had been assembled as
deliverable files. Seven of them (S4-S10) are already computed and sit in
``processed/``; S3 aggregates the per-cell entropy scores; S1 and S2 are curation
tables whose content lives in the manuscript's Methods and in ``msep.pathways``.

This script assembles all ten, writing one CSV per table plus a single workbook for
submission.

Usage::

    python notebooks/paper_figures/build_supplementary_tables.py

Writes to ``figures/Supplementary/tables/``:
    Table_S1.csv ... Table_S10.csv
    Supplementary_Tables_S1-S10.xlsx     (one sheet per table)

Two columns need author input before submission and are emitted empty:
S1 ``literature_source`` (per-gene citations; only pathway-level provenance is
recorded in Methods 2.4) and S2 ``reference``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

import msep.pathways as pathways

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
PROCESSED = REPO_ROOT / "processed"
DATA = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "figures" / "Supplementary" / "tables"

# Pathway-level provenance as stated in Methods 2.4. Per-gene citations are not
# recorded anywhere in the analysis output and must be supplied by the authors.
PATHWAY_PROVENANCE = {
    "ferroptosis": "FerrDb V2 and literature curation",
    "immune_evasion": "Literature curation",
    "emt": "Literature curation",
    "housekeeping": "Established reference gene panels",
}

# Functional sub-groupings, following the prose of Methods 2.4.
GENE_ROLES = {
    "ferroptosis": {
        "core regulator": ["GPX4", "ACSL4", "SLC7A11"],
        "iron metabolism": ["FTH1", "FTL", "NCOA4", "TFRC"],
        "lipid peroxidation": ["ALOX5", "ALOX12", "ALOX15", "FADS2", "SCD"],
        "antioxidant response": ["NFE2L2", "GCLC", "GCLM", "GSS"],
    },
    "immune_evasion": {
        "HLA class I": ["HLA-A", "HLA-B", "HLA-C"],
        "non-classical HLA": ["HLA-E", "HLA-F", "HLA-G"],
        "antigen presentation": ["B2M", "TAP1", "TAP2"],
        "NK ligand/receptor": ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
                               "KLRC1", "KLRD1"],
        "checkpoint ligand": ["CD274", "PDCD1LG2", "CD47", "LGALS9", "VTCN1"],
        "metalloproteinase": ["ADAM10", "ADAM17"],
    },
    "emt": {
        "mesenchymal marker": ["VIM", "CDH2", "FN1", "ACTA2"],
        "epithelial marker": ["CDH1", "KRT8", "KRT18", "KRT19"],
        "EMT transcription factor": ["ZEB1", "ZEB2", "SNAI1", "SNAI2",
                                     "TWIST1", "TWIST2"],
        "notochordal lineage": ["TBXT"],
        "TGF-beta component": ["TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
                               "SMAD2", "SMAD3", "SMAD4"],
    },
}

# Marker genes and criteria as given in Methods 2.3, with counts from Results 3.1.
CELL_TYPE_ANNOTATION = [
    ("T_cell", "CD3D, CD3E, CD8A", "Marker-based Leiden cluster assignment", 26644),
    ("Macrophage", "CD68, CD163, CSF1R", "Marker-based Leiden cluster assignment", 23203),
    ("NK_cell", "NKG7, GNLY, KLRD1", "Marker-based Leiden cluster assignment", 16931),
    ("Monocyte", "S100A8, S100A9, CD14", "Marker-based Leiden cluster assignment", 12232),
    ("Mast_cell", "KIT, TPSAB1", "Marker-based Leiden cluster assignment", 9054),
    ("CSC_TBXT+", "TBXT", "TBXT-positive malignant cells", 6730),
    ("B_cell", "CD79A, MS4A1", "Marker-based Leiden cluster assignment", 5462),
    ("Malignant_TBXT-", "TBXT (absent)", "TBXT-negative malignant cells", 2485),
    ("DC", "CLEC10A, CD1C", "Marker-based Leiden cluster assignment", 2431),
    ("Stromal", "DCN, COL1A1", "Marker-based Leiden cluster assignment", 1412),
]

ENTROPY_METRICS = [
    "entropy_global",
    "entropy_no_ribo",
    "entropy_ferroptosis_normalized",
    "entropy_immune_evasion_normalized",
    "entropy_emt_mechanical_normalized",
    "entropy_housekeeping_normalized",
]

TABLE_TITLES = {
    "S1": "Complete pathway gene sets with literature sources",
    "S2": "Cell type marker genes and annotation criteria",
    "S3": "Per-cell entropy statistics for all ten cell types",
    "S4": "Full pairwise Mann-Whitney U test results (42 comparisons, FDR-corrected)",
    "S5": "Pan-cancer across-cells CV (12 cancer types x 4 pathways)",
    "S6": "Gene-level CV for 23 key defense genes across 7 cell types",
    "S7": "Pseudo-perturbation full results (6 shield genes x 3 target pathways)",
    "S8": "XBP1 pan-cancer consolidation full statistics (9 cancer types x 3 pathways)",
    "S9": "Gene set sensitivity analysis, including overlap-only comparison",
    "S10": "Bootstrap 95% confidence intervals for pathway CV with overlap assessment",
}


def _read(name: str, base: Path = PROCESSED) -> pd.DataFrame:
    """Read a source CSV, with a pointed error if it is absent."""
    path = base / name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Supplementary tables are assembled from the full "
            "analysis output, which is not tracked in git."
        )
    return pd.read_csv(path)


def table_s1() -> pd.DataFrame:
    """Pathway gene sets, one row per gene."""
    gene_sets = {
        "ferroptosis": pathways.FERROPTOSIS,
        "immune_evasion": pathways.IMMUNE_EVASION,
        "emt": pathways.EMT,
        "housekeeping": pathways.HOUSEKEEPING,
    }
    rows = []
    for pathway, genes in gene_sets.items():
        roles = GENE_ROLES.get(pathway, {})
        gene_to_role = {g: role for role, gs in roles.items() for g in gs}
        for gene in genes:
            rows.append({
                "pathway": pathway,
                "gene": gene,
                "functional_role": gene_to_role.get(gene, ""),
                "set_provenance": PATHWAY_PROVENANCE[pathway],
                "literature_source": "",  # per-gene citation: author input required
            })
    return pd.DataFrame(rows)


def table_s2() -> pd.DataFrame:
    """Cell type annotation criteria."""
    df = pd.DataFrame(
        CELL_TYPE_ANNOTATION,
        columns=["cell_type", "marker_genes", "annotation_criterion", "n_cells"],
    )
    df["pct_of_total"] = (df["n_cells"] / df["n_cells"].sum() * 100).round(1)
    df["reference"] = ""  # author input required
    return df


def table_s3() -> pd.DataFrame:
    """Per-cell entropy statistics aggregated by cell type."""
    scores = _read("entropy_scores.csv")
    grouped = scores.groupby("cell_type_fine")

    out = pd.DataFrame({"n_cells": grouped.size()})
    for metric in ENTROPY_METRICS:
        if metric not in scores.columns:
            logger.warning("entropy_scores.csv has no column %s; skipping", metric)
            continue
        stats = grouped[metric]
        short = metric.replace("entropy_", "").replace("_normalized", "_norm")
        out[f"{short}_mean"] = stats.mean().round(4)
        out[f"{short}_median"] = stats.median().round(4)
        out[f"{short}_sd"] = stats.std().round(4)

    return out.sort_values("n_cells", ascending=False).reset_index()


def table_s9() -> pd.DataFrame:
    """Gene set sensitivity, merging the per-cell-type and overlap-only analyses."""
    extended = _read("gene_set_sensitivity_extended.csv")
    overlap_extended = _read("gene_set_sensitivity_overlap_extended.csv")
    summary = _read("gene_set_sensitivity.csv")
    overlap = _read("gene_set_sensitivity_overlap.csv")

    df = extended.merge(overlap_extended, on=["cell_type", "pathway"], how="outer")

    # Attach the pathway-level metadata: which alternative set was used, how big it
    # was, and how many genes the two definitions share.
    meta = summary[[
        "pathway", "original_source", "alternative_source",
        "original_n_total", "original_n_detected",
        "alternative_n_total", "alternative_n_detected", "overlap_n",
    ]].merge(overlap[["pathway", "overlap_genes"]], on="pathway", how="left")

    return df.merge(meta, on="pathway", how="left")


def build_all() -> dict[str, pd.DataFrame]:
    """Assemble every supplementary table."""
    return {
        "S1": table_s1(),
        "S2": table_s2(),
        "S3": table_s3(),
        "S4": _read("entropy_statistical_tests.csv"),
        "S5": _read("pan_cancer_cv.csv", base=DATA),
        "S6": _read("gene_level_cv_results.csv"),
        "S7": _read("perturbation_pseudo_results.csv"),
        "S8": _read("pan_cancer_xbp1_consolidation.csv"),
        "S9": table_s9(),
        "S10": _read("bootstrap_ci_overlap_table.csv"),
    }


def main() -> None:
    """Build all ten tables and write CSVs plus a combined workbook."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    tables = build_all()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(OUTPUT_DIR / f"Table_{name}.csv", index=False)

    workbook = OUTPUT_DIR / "Supplementary_Tables_S1-S10.xlsx"
    try:
        with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
            for name, df in tables.items():
                df.to_excel(writer, sheet_name=name, index=False)
        logger.info("Wrote %s", workbook.name)
    except ImportError:
        logger.warning("openpyxl not installed; CSVs written, workbook skipped")

    logger.info("Wrote %d tables to %s", len(tables), OUTPUT_DIR)
    for name, df in tables.items():
        logger.info("  Table %-4s %3d rows x %2d cols   %s",
                    name, len(df), df.shape[1], TABLE_TITLES[name])


if __name__ == "__main__":
    main()
