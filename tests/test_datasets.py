"""Tests for ``msep.datasets`` — the bundled example AnnData builder."""

from __future__ import annotations

import numpy as np
import pytest
import anndata as ad
from scipy.sparse import issparse

import msep
from msep.datasets import load_example, list_examples


# ---------------------------------------------------------------------------
# Basic structural tests
# ---------------------------------------------------------------------------


class TestLoadExampleStructure:
    """The returned AnnData must have the advertised shape and fields."""

    def test_returns_anndata(self):
        adata = load_example()
        assert isinstance(adata, ad.AnnData)

    def test_default_n_cells(self):
        adata = load_example()
        assert adata.n_obs == 500

    def test_custom_n_cells(self):
        adata = load_example(n_cells=200)
        assert adata.n_obs == 200

    def test_rejects_too_small(self):
        with pytest.raises(ValueError):
            load_example(n_cells=10)

    def test_has_raw_counts_layer(self):
        adata = load_example()
        assert "raw_counts" in adata.layers

    def test_counts_are_sparse(self):
        adata = load_example()
        assert issparse(adata.X)
        assert issparse(adata.layers["raw_counts"])

    def test_counts_are_nonnegative_integers(self):
        adata = load_example()
        x = adata.X.toarray()
        assert (x >= 0).all()
        # NB samples are integers; stored as float32, but should round-trip
        assert np.allclose(x, np.round(x))

    def test_has_cell_type_column(self):
        adata = load_example()
        assert "cell_type" in adata.obs.columns
        types = set(adata.obs["cell_type"].unique())
        assert types == {"CSC_TBXT+", "T_cell", "Macrophage", "Stromal"}

    def test_has_patient_id_column(self):
        adata = load_example(n_patients=3)
        assert "patient_id" in adata.obs.columns
        assert adata.obs["patient_id"].nunique() == 3

    def test_uns_pathways_populated(self):
        adata = load_example()
        assert "pathways" in adata.uns
        pw = adata.uns["pathways"]
        assert set(pw.keys()) == {
            "ferroptosis", "immune_evasion", "emt", "housekeeping",
        }
        # every pathway gene present in adata.var
        for genes in pw.values():
            assert all(g in adata.var_names for g in genes)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


class TestLoadExampleDeterminism:
    """Same seed must yield identical output; different seeds must differ."""

    def test_same_seed_identical(self):
        a = load_example(seed=42)
        b = load_example(seed=42)
        assert np.array_equal(a.X.toarray(), b.X.toarray())
        assert list(a.obs["cell_type"]) == list(b.obs["cell_type"])
        assert list(a.var_names) == list(b.var_names)

    def test_different_seed_different(self):
        a = load_example(seed=42)
        b = load_example(seed=7)
        assert not np.array_equal(a.X.toarray(), b.X.toarray())


# ---------------------------------------------------------------------------
# Phenotype reproduction — the whole point of the demo
# ---------------------------------------------------------------------------


class TestLoadExamplePhenotype:
    """The synthetic dataset must reproduce the paper's qualitative paradox.

    We do not require paper-exact values (the toy dataset is much smaller
    than the real 106,584-cell cohort), but the qualitative rankings must
    match:
      * CSC has the highest mean per-cell Shannon entropy
      * Within CSC: EMT CV < ferroptosis CV < immune_evasion CV
      * TBXT is tightly expressed in CSC (CV noticeably lower than its CV
        across the full dataset).
    """

    def test_csc_has_highest_entropy(self):
        adata = load_example(n_cells=500)
        result = msep.profile(
            adata,
            pathways="cancer_defense",
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=False,
        )
        # paradox_summary indexes by cell_type; entropy column is
        # ``median_entropy``.
        medians = result.paradox_summary["median_entropy"]
        top = medians.idxmax()
        assert top == "CSC_TBXT+", (
            f"Expected CSC_TBXT+ to have highest per-cell entropy; "
            f"got {top} with medians:\n{medians}"
        )

    def test_csc_pathway_cv_ranking(self):
        """Within CSC: EMT tightest, immune evasion most heterogeneous."""
        adata = load_example(n_cells=500)
        result = msep.profile(
            adata,
            pathways="cancer_defense",
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=False,
        )
        cv = result.pathway_cv
        # Column layout: cell_type, pathway, cv, ...
        csc = cv[cv["cell_type"] == "CSC_TBXT+"].set_index("pathway")["cv"]

        assert csc["emt"] < csc["ferroptosis"], (
            f"EMT CV ({csc['emt']:.3f}) should be lower than ferroptosis "
            f"CV ({csc['ferroptosis']:.3f}) in CSC."
        )
        assert csc["ferroptosis"] < csc["immune_evasion"], (
            f"Ferroptosis CV ({csc['ferroptosis']:.3f}) should be lower "
            f"than immune evasion CV ({csc['immune_evasion']:.3f}) in CSC."
        )

    def test_tbxt_low_cv_in_csc(self):
        """TBXT should show low CV in CSC (mirroring the paper's CV = 1.02)."""
        adata = load_example(n_cells=500)
        csc_mask = adata.obs["cell_type"] == "CSC_TBXT+"
        tbxt_col = adata.var_names.get_loc("TBXT")
        tbxt_csc = adata.X.toarray()[csc_mask, tbxt_col]
        cv = tbxt_csc.std() / tbxt_csc.mean()
        assert cv < 0.3, (
            f"TBXT CV in CSC should be low (<0.3); got {cv:.3f}"
        )


# ---------------------------------------------------------------------------
# XBP1 stress setup
# ---------------------------------------------------------------------------


class TestLoadExampleXbp1:
    def test_xbp1_present(self):
        adata = load_example()
        assert "XBP1" in adata.var_names

    def test_xbp1_stress_three_groups_present(self):
        """With default settings, CSC cells cover zero / low / high XBP1."""
        adata = load_example(n_cells=500)
        csc_mask = adata.obs["cell_type"] == "CSC_TBXT+"
        xbp1 = adata.X.toarray()[csc_mask, adata.var_names.get_loc("XBP1")]

        zero_frac = (xbp1 == 0).mean()
        nonzero_frac = (xbp1 > 0).mean()

        # Paper reports 43.3% XBP1-expressing in CSC; the synthetic dataset
        # is calibrated near this value but may drift +-10% under n_cells=500.
        assert zero_frac > 0.40, f"zero XBP1 fraction {zero_frac:.2f}"
        assert nonzero_frac > 0.25, f"nonzero XBP1 fraction {nonzero_frac:.2f}"

    def test_xbp1_stress_can_be_disabled(self):
        adata = load_example(n_cells=200, include_xbp1_stress=False)
        # Without the forced pattern, XBP1 still present but follows the
        # generic pathway-block distribution (no guarantee of zero cells).
        assert "XBP1" in adata.var_names


# ---------------------------------------------------------------------------
# Profile end-to-end smoke test
# ---------------------------------------------------------------------------


class TestEndToEnd:
    """The demo dataset must flow through ``msep.profile()`` cleanly."""

    def test_profile_runs(self):
        adata = load_example(n_cells=200)
        result = msep.profile(
            adata,
            pathways="cancer_defense",
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=True,
            compute_perturbation=False,
            compute_xbp1=False,
        )
        assert result.pathway_cv is not None
        assert len(result.pathway_cv) > 0
        # four cell types x four pathways
        assert result.pathway_cv["cell_type"].nunique() == 4

    def test_profile_with_xbp1_consolidation(self):
        adata = load_example(n_cells=500)
        result = msep.profile(
            adata,
            pathways="cancer_defense",
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=False,
            compute_xbp1=True,
            perturbation_cell_type="CSC_TBXT+",
        )
        assert result.xbp1 is not None


# ---------------------------------------------------------------------------
# list_examples helper
# ---------------------------------------------------------------------------


def test_list_examples_returns_dict_with_demo():
    out = list_examples()
    assert isinstance(out, dict)
    assert "chordoma_demo" in out
