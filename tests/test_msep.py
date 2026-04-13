"""
Tests for the msep package.

Uses synthetic AnnData objects so tests run without real data.
"""

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix

import msep
from msep.entropy import per_cell_entropy, normalized_entropy
from msep.coordination import pathway_cv, pathway_cv_table, bootstrap_cv, fano_factor, gene_level_cv
from msep.perturbation import stratify_by_gene, pseudo_perturbation, xbp1_consolidation
from msep.pathways import get_pathways, resolve_genes, CANCER_DEFENSE


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def synthetic_adata():
    """Create a small synthetic AnnData with known properties."""
    np.random.seed(42)

    n_cells = 500
    n_genes = 100

    gene_names = [f"GENE{i}" for i in range(n_genes)]
    # Add some known pathway genes
    pathway_genes = ["VIM", "CDH1", "FN1", "TBXT", "SMAD2", "SMAD3",
                     "GPX4", "SLC7A11", "FTH1", "ACSL4",
                     "HLA-E", "B2M", "CD274", "MICA",
                     "ACTB", "GAPDH", "XBP1"]
    gene_names[:len(pathway_genes)] = pathway_genes

    # Generate count data: Poisson with varying rates
    rates = np.random.exponential(2.0, size=(1, n_genes))
    counts = np.random.poisson(rates, size=(n_cells, n_genes)).astype(np.float32)

    # Make VIM uniformly expressed (low CV) in first 200 cells
    counts[:200, gene_names.index("VIM")] = np.random.poisson(10, 200)
    # Make HLA-E variable (high CV)
    counts[:200, gene_names.index("HLA-E")] = np.random.poisson(
        np.random.exponential(5, 200)
    )

    # XBP1: zero in first 250, low in next 150, high in last 100
    counts[:250, gene_names.index("XBP1")] = 0
    counts[250:400, gene_names.index("XBP1")] = np.random.poisson(2, 150)
    counts[400:, gene_names.index("XBP1")] = np.random.poisson(15, 100)

    cell_types = np.array(
        ["CSC"] * 200 + ["T_cell"] * 150 + ["Macrophage"] * 150
    )

    adata = ad.AnnData(
        X=csr_matrix(counts),
        obs=pd.DataFrame({"cell_type": cell_types},
                         index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=gene_names),
    )
    adata.layers["raw_counts"] = adata.X.copy()
    return adata


@pytest.fixture
def small_pathways():
    """Minimal pathway dict for testing."""
    return {
        "emt": ["VIM", "CDH1", "FN1", "TBXT", "SMAD2", "SMAD3"],
        "ferroptosis": ["GPX4", "SLC7A11", "FTH1", "ACSL4"],
        "immune_evasion": ["HLA-E", "B2M", "CD274", "MICA"],
        "housekeeping": ["ACTB", "GAPDH"],
    }


# -----------------------------------------------------------------------
# Pathway tests
# -----------------------------------------------------------------------

class TestPathways:
    def test_cancer_defense_has_four_sets(self):
        pw = get_pathways("cancer_defense")
        assert set(pw.keys()) == {"ferroptosis", "immune_evasion", "emt", "housekeeping"}

    def test_custom_pathways(self):
        custom = {"my_set": ["GENE1", "GENE2"]}
        pw = get_pathways(custom=custom)
        assert pw == custom

    def test_invalid_custom_raises(self):
        with pytest.raises(ValueError):
            get_pathways(custom={"empty": []})

    def test_resolve_genes(self, synthetic_adata):
        present, missing, idx = resolve_genes(
            synthetic_adata.var_names, ["VIM", "CDH1", "NONEXISTENT"]
        )
        assert "VIM" in present
        assert "NONEXISTENT" in missing
        assert len(idx) == 2


# -----------------------------------------------------------------------
# Entropy tests
# -----------------------------------------------------------------------

class TestEntropy:
    def test_shape(self, synthetic_adata):
        mat = synthetic_adata.layers["raw_counts"]
        h, n = per_cell_entropy(mat)
        assert h.shape == (500,)
        assert n.shape == (500,)

    def test_non_negative(self, synthetic_adata):
        mat = synthetic_adata.layers["raw_counts"]
        h, _ = per_cell_entropy(mat)
        valid = ~np.isnan(h)
        assert np.all(h[valid] >= 0)

    def test_max_entropy_bound(self, synthetic_adata):
        """Entropy cannot exceed log2(n_expressed)."""
        mat = synthetic_adata.layers["raw_counts"]
        h, n = per_cell_entropy(mat)
        valid = (~np.isnan(h)) & (n > 1)
        max_possible = np.log2(n[valid].astype(float))
        assert np.all(h[valid] <= max_possible + 1e-10)

    def test_uniform_distribution_max(self):
        """Uniform counts → maximum entropy."""
        mat = np.ones((1, 100), dtype=float) * 10
        h, n = per_cell_entropy(mat)
        expected = np.log2(100)
        np.testing.assert_allclose(h[0], expected, atol=1e-10)

    def test_single_gene_zero_entropy(self):
        """All reads in one gene → entropy = 0."""
        mat = np.zeros((1, 100), dtype=float)
        mat[0, 0] = 1000
        h, _ = per_cell_entropy(mat)
        assert h[0] == 0.0

    def test_sparse_dense_equal(self, synthetic_adata):
        """Sparse and dense inputs give same result."""
        mat_sparse = synthetic_adata.layers["raw_counts"]
        mat_dense = mat_sparse.toarray()
        h_sparse, _ = per_cell_entropy(mat_sparse)
        h_dense, _ = per_cell_entropy(mat_dense)
        valid = ~np.isnan(h_sparse)
        np.testing.assert_allclose(h_sparse[valid], h_dense[valid], atol=1e-10)

    def test_normalized_entropy_range(self, synthetic_adata):
        mat = synthetic_adata.layers["raw_counts"]
        h, n = per_cell_entropy(mat)
        h_norm = normalized_entropy(h, n)
        valid = ~np.isnan(h_norm)
        assert np.all(h_norm[valid] >= 0)
        assert np.all(h_norm[valid] <= 1 + 1e-10)

    def test_gene_subset(self, synthetic_adata):
        mat = synthetic_adata.layers["raw_counts"]
        idx = np.array([0, 1, 2])
        h_sub, _ = per_cell_entropy(mat, gene_indices=idx)
        h_full, _ = per_cell_entropy(mat)
        # Subset entropy ≤ full entropy
        valid = ~np.isnan(h_sub) & ~np.isnan(h_full)
        assert np.all(h_sub[valid] <= h_full[valid] + 1e-10)


# -----------------------------------------------------------------------
# Coordination tests
# -----------------------------------------------------------------------

class TestCoordination:
    def test_pathway_cv_basic(self, synthetic_adata, small_pathways):
        mat = synthetic_adata.layers["raw_counts"]
        cv, ng = pathway_cv(
            mat, synthetic_adata.var_names, small_pathways["emt"]
        )
        assert not np.isnan(cv)
        assert cv > 0
        assert ng == 6

    def test_pathway_cv_missing_genes(self, synthetic_adata):
        cv, ng = pathway_cv(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            ["NONEXISTENT1", "NONEXISTENT2"],
        )
        assert np.isnan(cv)
        assert ng == 0

    def test_cv_table(self, synthetic_adata, small_pathways):
        df = pathway_cv_table(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            small_pathways,
            cell_type_labels=np.array(synthetic_adata.obs["cell_type"]),
        )
        assert "cell_type" in df.columns
        assert "pathway" in df.columns
        assert len(df) == 3 * 4  # 3 cell types × 4 pathways

    def test_bootstrap(self, synthetic_adata, small_pathways):
        result = bootstrap_cv(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            small_pathways["emt"],
            n_boot=100,
        )
        assert result["ci_low"] < result["ci_high"]
        assert not np.isnan(result["cv_mean"])

    def test_fano(self, synthetic_adata, small_pathways):
        f, ng = fano_factor(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            small_pathways["emt"],
        )
        assert not np.isnan(f)
        assert f > 0

    def test_gene_level_cv(self, synthetic_adata, small_pathways):
        df = gene_level_cv(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            small_pathways["emt"],
        )
        assert "gene" in df.columns
        assert "cv" in df.columns
        assert len(df) == 6

    def test_uniform_low_cv(self):
        """Uniformly expressed gene → CV = 0."""
        mat = np.full((100, 3), 10.0)
        var_names = pd.Index(["A", "B", "C"])
        cv, _ = pathway_cv(mat, var_names, ["A", "B", "C"])
        assert cv == 0.0


# -----------------------------------------------------------------------
# Perturbation tests
# -----------------------------------------------------------------------

class TestPerturbation:
    def test_stratify(self, synthetic_adata):
        strat = stratify_by_gene(
            synthetic_adata.layers["raw_counts"],
            synthetic_adata.var_names,
            "VIM",
        )
        assert strat["n_high"] > 0
        assert strat["n_low"] > 0
        # no overlap
        assert not np.any(strat["high_mask"] & strat["low_mask"])

    def test_pseudo_perturbation(self, synthetic_adata, small_pathways):
        mat = synthetic_adata.layers["raw_counts"]
        df = pseudo_perturbation(
            mat, synthetic_adata.var_names,
            shield_genes=["VIM", "GPX4"],
            pathways=small_pathways,
            n_perm=50,  # small for speed
        )
        assert len(df) > 0
        assert "delta_cv" in df.columns
        assert "perm_p" in df.columns

    def test_xbp1(self, synthetic_adata, small_pathways):
        mat = synthetic_adata.layers["raw_counts"]
        df = xbp1_consolidation(mat, synthetic_adata.var_names, small_pathways)
        assert len(df) > 0
        groups = set(df["xbp1_group"])
        assert "XBP1-zero" in groups
        assert "XBP1-high" in groups


# -----------------------------------------------------------------------
# Integration test: profile()
# -----------------------------------------------------------------------

class TestProfile:
    def test_basic_profile(self, synthetic_adata, small_pathways):
        result = msep.profile(
            synthetic_adata,
            pathways=small_pathways,
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=True,
        )
        assert isinstance(result, msep.MSEPResult)
        assert len(result.per_cell_entropy) == 500
        assert len(result.pathway_cv) == 12  # 3 types × 4 pathways

    def test_paradox_summary(self, synthetic_adata, small_pathways):
        result = msep.profile(
            synthetic_adata,
            pathways=small_pathways,
            cell_type_key="cell_type",
            compute_bootstrap=False,
        )
        summary = result.paradox_summary
        assert "median_entropy" in summary.columns
        assert len(summary) == 3

    def test_full_pipeline(self, synthetic_adata, small_pathways):
        result = msep.profile(
            synthetic_adata,
            pathways=small_pathways,
            cell_type_key="cell_type",
            compute_bootstrap=True,
            compute_gene_cv=True,
            compute_perturbation=True,
            compute_xbp1=True,
            shield_genes=["VIM", "GPX4"],
            n_boot=50,
            n_perm=20,
        )
        assert result.bootstrap is not None
        assert result.perturbation is not None
        assert result.xbp1 is not None
        assert result.gene_cv is not None

    def test_repr(self, synthetic_adata, small_pathways):
        result = msep.profile(
            synthetic_adata,
            pathways=small_pathways,
            cell_type_key="cell_type",
            compute_bootstrap=False,
        )
        r = repr(result)
        assert "MSEPResult" in r
        assert "500" in r

    def test_builtin_pathways(self, synthetic_adata):
        """profile() works with built-in 'cancer_defense' pathways."""
        result = msep.profile(
            synthetic_adata,
            pathways="cancer_defense",
            cell_type_key="cell_type",
            compute_bootstrap=False,
            compute_gene_cv=False,
        )
        assert len(result.pathway_cv) > 0
