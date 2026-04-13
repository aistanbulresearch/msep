"""
Built-in pathway gene sets for multi-scale entropy profiling.

Gene sets curated from FerrDb V2, KEGG, Reactome, and literature.
Users may supply their own gene sets as ``dict[str, list[str]]``.
"""

from __future__ import annotations

from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# Curated cancer defense pathway gene sets
# ---------------------------------------------------------------------------

FERROPTOSIS: List[str] = [
    "GPX4", "ACSL4", "SLC7A11", "FTH1", "FTL", "NCOA4", "LPCAT3",
    "TFRC", "GSS", "GCLC", "NFE2L2", "HMOX1", "SCD", "FADS2", "CBS",
    "GCLM", "SLC3A2", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",
    "ACSL3", "CARS1", "CISD1", "CISD2", "DPP4", "FANCD2", "GLS2",
    "MT1G", "PEBP1", "SAT1", "VDAC2", "VDAC3",
]

IMMUNE_EVASION: List[str] = [
    "HLA-E", "B2M", "CD274", "MICA", "MICB", "ADAM10", "ADAM17",
    "TAP1", "TAP2", "KLRC1", "KLRD1", "ULBP1", "ULBP2", "ULBP3",
    "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLA-G",
    "CD47", "LGALS9", "PVRL2", "PVR", "VTCN1", "TNFRSF14",
    "NECTIN2", "CD200", "IDO1", "PDCD1LG2",
]

EMT: List[str] = [
    "VIM", "CDH1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2",
    "TWIST1", "TWIST2", "KRT18", "KRT8", "KRT19", "TBXT", "FN1",
    "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
    "SMAD2", "SMAD3", "SMAD4", "ACTA2",
    "MMP2", "MMP9", "SPARC", "ITGB1", "ITGAV",
]

HOUSEKEEPING: List[str] = [
    "ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "RPLP0",
    "PGK1", "TFRC", "B2M", "GUSB", "HMBS", "YWHAZ",
    "SDHA", "UBC", "HSP90AB1", "LDHA", "ALDOA", "ENO1",
    "PKM", "TUBB",
]

# ---------------------------------------------------------------------------
# Convenience collections
# ---------------------------------------------------------------------------

CANCER_DEFENSE: Dict[str, List[str]] = {
    "ferroptosis": FERROPTOSIS,
    "immune_evasion": IMMUNE_EVASION,
    "emt": EMT,
    "housekeeping": HOUSEKEEPING,
}
"""Default pathway collection for cancer defense profiling."""


def get_pathways(
    name: str = "cancer_defense",
    custom: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, List[str]]:
    """Return a pathway dictionary by name, or validate a custom one.

    Parameters
    ----------
    name
        Built-in collection name.  Currently only ``"cancer_defense"``.
    custom
        If provided, returned as-is after validation.  Each value must
        be a non-empty list of gene-name strings.

    Returns
    -------
    dict
        Mapping of pathway name to gene list.

    Raises
    ------
    ValueError
        If *name* is unknown and *custom* is ``None``, or if any value
        in *custom* is empty or not a list of strings.
    """
    if custom is not None:
        _validate(custom)
        return custom

    collections = {
        "cancer_defense": CANCER_DEFENSE,
    }
    if name not in collections:
        raise ValueError(
            f"Unknown pathway collection {name!r}. "
            f"Available: {list(collections)}. "
            f"Pass custom=dict(...) for user-defined gene sets."
        )
    return collections[name]


def _validate(pathways: Dict[str, List[str]]) -> None:
    if not isinstance(pathways, dict) or not pathways:
        raise ValueError("pathways must be a non-empty dict")
    for key, genes in pathways.items():
        if not isinstance(genes, (list, tuple)) or len(genes) == 0:
            raise ValueError(
                f"Pathway {key!r} must be a non-empty list of gene names"
            )
        if not all(isinstance(g, str) for g in genes):
            raise ValueError(
                f"All entries in pathway {key!r} must be strings"
            )


def resolve_genes(
    var_names,
    gene_list: List[str],
) -> tuple:
    """Match a gene list against the dataset's var_names.

    Parameters
    ----------
    var_names
        ``adata.var_names`` (Index or array-like).
    gene_list
        Gene symbols to look up.

    Returns
    -------
    present : list[str]
        Genes found in *var_names*.
    missing : list[str]
        Genes not found.
    indices : ndarray
        Integer positions of *present* genes in *var_names*.
    """
    import numpy as np
    import pandas as pd

    if not isinstance(var_names, pd.Index):
        var_names = pd.Index(var_names)

    present = [g for g in gene_list if g in var_names]
    missing = [g for g in gene_list if g not in var_names]
    indices = np.where(var_names.isin(present))[0]
    return present, missing, indices
