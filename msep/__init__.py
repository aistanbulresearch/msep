"""
msep — Multi-Scale Entropy Profiling for single-cell transcriptomics
====================================================================

Integrates per-cell Shannon entropy with across-cells coefficient of
variation (CV), decomposed by biological pathway, to characterise
population-level transcriptomic coordination in cancer and beyond.

Quick start
-----------
>>> import msep
>>> result = msep.profile(adata, cell_type_key="cell_type")
>>> result.paradox_summary          # entropy vs CV per cell type
>>> result.pathway_cv               # pathway-level CV table

Reference
---------
Çavuş & Kuşkucu (2026).  Multi-Scale Entropy Profiling Reveals
Pathway-Selective Defense Coordination Across Cancer Types.
"""

__version__ = "1.1.0"

# Core API
from .core import profile  # noqa: F401
from .result import MSEPResult  # noqa: F401

# Sub-module access
from . import entropy  # noqa: F401
from . import coordination  # noqa: F401
from . import perturbation  # noqa: F401
from . import pathways  # noqa: F401
from . import datasets  # noqa: F401

# Plotting (lazy — only imported when called)
def plot_entropy_violin(result, **kw):
    """Violin plot of per-cell entropy.  See :mod:`msep.plotting`."""
    from .plotting import plot_entropy_violin as _fn
    return _fn(result, **kw)

def plot_pathway_cv_heatmap(result, **kw):
    """CV heatmap.  See :mod:`msep.plotting`."""
    from .plotting import plot_pathway_cv_heatmap as _fn
    return _fn(result, **kw)

def plot_paradox(result, **kw):
    """Paradox scatter.  See :mod:`msep.plotting`."""
    from .plotting import plot_paradox as _fn
    return _fn(result, **kw)

def plot_pan_cancer(cv_data, **kw):
    """Pan-cancer ranking.  See :mod:`msep.plotting`."""
    from .plotting import plot_pan_cancer as _fn
    return _fn(cv_data, **kw)

# Bayesian (lazy — requires scvi-tools)
def bayesian_validate(adata, pathways, **kw):
    """Bayesian variance decomposition via scVI.  See :mod:`msep.bayesian`."""
    from .bayesian import bayesian_variance_decomposition as _fn
    return _fn(adata, pathways, **kw)


__all__ = [
    "profile",
    "MSEPResult",
    "entropy",
    "coordination",
    "perturbation",
    "pathways",
    "datasets",
    "bayesian_validate",
    "plot_entropy_violin",
    "plot_pathway_cv_heatmap",
    "plot_paradox",
    "plot_pan_cancer",
]
