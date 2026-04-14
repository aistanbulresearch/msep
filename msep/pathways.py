"""
Built-in and external pathway gene sets for multi-scale entropy profiling.

Includes curated cancer defense gene sets and seamless access to MSigDB
(Hallmark, C2/KEGG, C2/Reactome, C5/GO, C6/Oncogenic), enabling
application across any biological domain.

Users may also supply their own gene sets as ``dict[str, list[str]]``.

Examples
--------
>>> import msep
>>> # Built-in curated sets
>>> result = msep.profile(adata, pathways="cancer_defense")
>>>
>>> # MSigDB Hallmark (50 gene sets, auto-downloaded)
>>> result = msep.profile(adata, pathways="hallmark")
>>>
>>> # KEGG pathways
>>> result = msep.profile(adata, pathways="kegg")
>>>
>>> # Pick specific MSigDB sets
>>> pw = msep.pathways.from_msigdb("hallmark",
...     include=["EPITHELIAL_MESENCHYMAL_TRANSITION", "INFLAMMATORY_RESPONSE"])
>>> result = msep.profile(adata, pathways=pw)
>>>
>>> # Combine MSigDB with custom genes
>>> pw = msep.pathways.from_msigdb("hallmark", top_n=5)
>>> pw["my_custom_set"] = ["GENE1", "GENE2", "GENE3"]
>>> result = msep.profile(adata, pathways=pw)
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union

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
# Built-in collections
# ---------------------------------------------------------------------------

CANCER_DEFENSE: Dict[str, List[str]] = {
    "ferroptosis": FERROPTOSIS,
    "immune_evasion": IMMUNE_EVASION,
    "emt": EMT,
    "housekeeping": HOUSEKEEPING,
}
"""Default pathway collection for cancer defense profiling."""

# MSigDB collection aliases
_MSIGDB_ALIASES = {
    "hallmark": "h.all",
    "kegg": "c2.cp.kegg_legacy",
    "reactome": "c2.cp.reactome",
    "biocarta": "c2.cp.biocarta",
    "pid": "c2.cp.pid",
    "wikipathways": "c2.cp.wikipathways",
    "go_bp": "c5.go.bp",
    "go_cc": "c5.go.cc",
    "go_mf": "c5.go.mf",
    "oncogenic": "c6.all",
    "immunologic": "c7.all",
    "cell_type": "c8.all",
}

# Cache directory
_CACHE_DIR = Path.home() / ".cache" / "msep" / "genesets"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def get_pathways(
    name: str = "cancer_defense",
    custom: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, List[str]]:
    """Return a pathway dictionary by name, or validate a custom one.

    Parameters
    ----------
    name
        Built-in name (``"cancer_defense"``) or MSigDB collection name
        (``"hallmark"``, ``"kegg"``, ``"reactome"``, ``"go_bp"``,
        ``"oncogenic"``, ``"immunologic"``, etc.).
        Use ``"hallmark:SET_NAME"`` to select a single set.
    custom
        If provided, returned as-is after validation.

    Returns
    -------
    dict
        Mapping of pathway name to gene list.
    """
    if custom is not None:
        _validate(custom)
        return custom

    # Built-in curated collections
    if name == "cancer_defense":
        return CANCER_DEFENSE

    # Single set from a collection: "hallmark:EPITHELIAL_MESENCHYMAL_TRANSITION"
    if ":" in name:
        collection, set_name = name.split(":", 1)
        all_sets = from_msigdb(collection)
        # Try exact match, then case-insensitive partial match
        if set_name in all_sets:
            return {set_name: all_sets[set_name]}
        matches = {k: v for k, v in all_sets.items()
                   if set_name.upper() in k.upper()}
        if matches:
            return matches
        raise ValueError(
            f"Gene set {set_name!r} not found in {collection!r}. "
            f"Available: {list(all_sets.keys())[:10]}..."
        )

    # MSigDB collections
    if name.lower() in _MSIGDB_ALIASES or name.startswith("c") or name == "h.all":
        return from_msigdb(name)

    raise ValueError(
        f"Unknown pathway collection {name!r}. "
        f"Available built-in: 'cancer_defense'. "
        f"MSigDB collections: {list(_MSIGDB_ALIASES.keys())}. "
        f"Or pass custom=dict(...)."
    )


# ---------------------------------------------------------------------------
# MSigDB integration
# ---------------------------------------------------------------------------

def from_msigdb(
    collection: str = "hallmark",
    species: str = "human",
    include: Optional[List[str]] = None,
    exclude: Optional[List[str]] = None,
    top_n: Optional[int] = None,
    min_genes: int = 5,
    max_genes: int = 500,
) -> Dict[str, List[str]]:
    """Fetch gene sets from MSigDB.

    Downloads GMT files from MSigDB and caches them locally.
    Requires internet connection on first use; cached thereafter.

    Parameters
    ----------
    collection
        Collection name. Shortcuts: ``"hallmark"`` (50 sets),
        ``"kegg"`` (186 sets), ``"reactome"`` (1,615 sets),
        ``"go_bp"`` (GO Biological Process), ``"go_cc"`` (GO Cellular
        Component), ``"go_mf"`` (GO Molecular Function),
        ``"oncogenic"`` (189 sets), ``"immunologic"`` (5,219 sets),
        ``"cell_type"`` (830 sets), ``"biocarta"``, ``"pid"``,
        ``"wikipathways"``.
    species
        ``"human"`` or ``"mouse"``.
    include
        If given, only return sets whose names contain any of these
        substrings (case-insensitive).
    exclude
        If given, exclude sets whose names contain any of these
        substrings (case-insensitive).
    top_n
        If given, return only the first *top_n* sets (alphabetically).
    min_genes, max_genes
        Filter sets by gene count.

    Returns
    -------
    dict
        ``{set_name: [gene_symbol, ...]}``.

    Examples
    --------
    >>> pw = msep.pathways.from_msigdb("hallmark")
    >>> len(pw)  # 50 Hallmark gene sets
    50

    >>> pw = msep.pathways.from_msigdb("kegg", include=["APOPTOSIS", "P53"])
    >>> list(pw.keys())
    ['KEGG_APOPTOSIS', 'KEGG_P53_SIGNALING_PATHWAY']

    >>> pw = msep.pathways.from_msigdb("reactome", max_genes=50, top_n=10)
    """
    # Resolve alias
    msigdb_name = _MSIGDB_ALIASES.get(collection.lower(), collection)

    # Try to load from cache first
    cache_file = _CACHE_DIR / f"{msigdb_name}.{species}.gmt"
    if cache_file.exists():
        gene_sets = _parse_gmt(cache_file)
    else:
        gene_sets = _download_msigdb(msigdb_name, species, cache_file)

    # Filter by gene count
    gene_sets = {
        k: v for k, v in gene_sets.items()
        if min_genes <= len(v) <= max_genes
    }

    # Include filter
    if include is not None:
        include_upper = [s.upper() for s in include]
        gene_sets = {
            k: v for k, v in gene_sets.items()
            if any(inc in k.upper() for inc in include_upper)
        }

    # Exclude filter
    if exclude is not None:
        exclude_upper = [s.upper() for s in exclude]
        gene_sets = {
            k: v for k, v in gene_sets.items()
            if not any(exc in k.upper() for exc in exclude_upper)
        }

    # Sort alphabetically and limit
    gene_sets = dict(sorted(gene_sets.items()))
    if top_n is not None:
        gene_sets = dict(list(gene_sets.items())[:top_n])

    if not gene_sets:
        warnings.warn(
            f"No gene sets found for collection={collection!r} with "
            f"current filters. Try adjusting include/exclude/min_genes/max_genes.",
            stacklevel=2,
        )

    return gene_sets


def list_collections() -> Dict[str, str]:
    """List available MSigDB collection shortcuts.

    Returns
    -------
    dict
        ``{shortcut: msigdb_identifier}``.
    """
    return dict(_MSIGDB_ALIASES)


def search_msigdb(
    query: str,
    collection: str = "hallmark",
    species: str = "human",
) -> Dict[str, List[str]]:
    """Search gene set names within a collection.

    Parameters
    ----------
    query
        Search term (case-insensitive substring match).
    collection
        MSigDB collection to search.
    species
        ``"human"`` or ``"mouse"``.

    Returns
    -------
    dict
        Matching gene sets.
    """
    all_sets = from_msigdb(collection, species=species)
    query_upper = query.upper()
    return {k: v for k, v in all_sets.items() if query_upper in k.upper()}


# ---------------------------------------------------------------------------
# GMT file handling
# ---------------------------------------------------------------------------

def _download_msigdb(
    msigdb_name: str,
    species: str,
    cache_file: Path,
) -> Dict[str, List[str]]:
    """Download GMT file from MSigDB."""
    import urllib.request
    import urllib.error

    species_tag = "Hs" if species == "human" else "Mm"
    url = (
        f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.{species_tag}/"
        f"{msigdb_name}.v2024.1.{species_tag}.symbols.gmt"
    )

    cache_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        urllib.request.urlretrieve(url, cache_file)
    except urllib.error.URLError:
        # Try alternative URL format
        url_alt = (
            f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.{species_tag}/"
            f"{msigdb_name}.v2023.2.{species_tag}.symbols.gmt"
        )
        try:
            urllib.request.urlretrieve(url_alt, cache_file)
        except urllib.error.URLError as e:
            raise ConnectionError(
                f"Could not download MSigDB collection {msigdb_name!r}. "
                f"Check your internet connection or collection name. "
                f"Available shortcuts: {list(_MSIGDB_ALIASES.keys())}. "
                f"Error: {e}"
            ) from e

    return _parse_gmt(cache_file)


def _parse_gmt(path: Path) -> Dict[str, List[str]]:
    """Parse a GMT (Gene Matrix Transposed) file.

    GMT format: each line is tab-separated with:
    SET_NAME \\t DESCRIPTION \\t GENE1 \\t GENE2 \\t ...
    """
    gene_sets = {}
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            if genes:
                gene_sets[name] = genes
    return gene_sets


def load_gmt(path: Union[str, Path]) -> Dict[str, List[str]]:
    """Load gene sets from a local GMT file.

    Parameters
    ----------
    path
        Path to a ``.gmt`` file.

    Returns
    -------
    dict
        ``{set_name: [gene_symbol, ...]}``.
    """
    return _parse_gmt(Path(path))


# ---------------------------------------------------------------------------
# Validation and gene resolution
# ---------------------------------------------------------------------------

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
