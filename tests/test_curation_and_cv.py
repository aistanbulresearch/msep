"""Regression tests for two curation/computation fixes.

1. ``coordination.pathway_cv`` averaged CV over genes with a positive mean but
   reported ``n_genes`` as the count of genes merely *found*, so the reported gene
   count could exceed the number of genes actually contributing to the CV.

2. ``pathways.IMMUNE_EVASION`` listed both ``PVRL2`` and ``NECTIN2``, which are two
   symbols for the same gene (PVRL2 is the withdrawn HGNC alias of NECTIN2). The set
   was declared as 29 genes but only ever resolved to 28 distinct genes.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import msep.pathways as pathways
from msep.coordination import pathway_cv

DATA_CSV = Path(__file__).resolve().parents[1] / "data" / "pan_cancer_cv.csv"


class TestPathwayCvGeneCount:
    """n_genes must report the genes actually averaged, not merely found."""

    def test_reports_genes_actually_averaged(self) -> None:
        # Four pathway genes are present in the matrix, but one is zero in every
        # cell (mean 0), so its CV is undefined and it is dropped from the average.
        # n_genes must therefore be 3, not 4.
        genes = ["G1", "G2", "G3", "G4"]
        var_names = np.array(genes)
        rng = np.random.default_rng(0)
        counts = rng.poisson(5, size=(200, 4)).astype(float)
        counts[:, 3] = 0.0  # G4 never expressed

        cv, n_genes = pathway_cv(counts, var_names, genes)

        assert np.isfinite(cv)
        assert n_genes == 3

    def test_all_expressed_reports_full_count(self) -> None:
        genes = ["G1", "G2", "G3"]
        var_names = np.array(genes)
        rng = np.random.default_rng(1)
        counts = rng.poisson(5, size=(150, 3)).astype(float) + 1.0

        _, n_genes = pathway_cv(counts, var_names, genes)
        assert n_genes == 3

    def test_absent_genes_not_counted(self) -> None:
        # Only two of the three requested genes exist in the matrix.
        var_names = np.array(["G1", "G2"])
        rng = np.random.default_rng(2)
        counts = rng.poisson(5, size=(100, 2)).astype(float) + 1.0

        _, n_genes = pathway_cv(counts, var_names, ["G1", "G2", "G3"])
        assert n_genes == 2


class TestImmuneEvasionGeneSet:
    """PVRL2 and NECTIN2 are the same gene; only one belongs in the set."""

    def test_no_pvrl2_nectin2_duplication(self) -> None:
        assert not ("PVRL2" in pathways.IMMUNE_EVASION
                    and "NECTIN2" in pathways.IMMUNE_EVASION)

    def test_nectin2_is_the_retained_symbol(self) -> None:
        # NECTIN2 is the current approved HGNC symbol; keep it, drop the alias.
        assert "NECTIN2" in pathways.IMMUNE_EVASION
        assert "PVRL2" not in pathways.IMMUNE_EVASION

    def test_immune_evasion_has_28_distinct_genes(self) -> None:
        assert len(pathways.IMMUNE_EVASION) == 28
        assert len(set(pathways.IMMUNE_EVASION)) == 28

    def test_other_pathways_unchanged(self) -> None:
        assert len(pathways.FERROPTOSIS) == 33
        assert len(pathways.EMT) == 28
        assert len(pathways.HOUSEKEEPING) == 20


class TestFerroptosisGeneSet:
    """The ferroptosis set must match the one the manuscript's figures were built on.

    Two different 33-gene ferroptosis sets were in circulation: the one the package
    originally shipped, and the one the Colab analysis (and therefore Figure 2, the
    scVI decomposition, and gene-set-sensitivity) actually used. They shared only 28
    genes. The analysis set is authoritative -- it is what the published figures
    report -- so the package is standardized to it here, and Supplementary Table S1
    and the pan-cancer ferroptosis values follow.
    """

    # The five genes the analysis set includes that the old package set omitted.
    ANALYSIS_ONLY = {"HSPB1", "IREB2", "NFS1", "PROM2", "STEAP3"}
    # The five the old package set carried that the analysis set does not.
    REMOVED = {"ALOX15B", "CARS1", "CISD2", "PEBP1", "SLC3A2"}

    def test_has_33_distinct_genes(self) -> None:
        assert len(pathways.FERROPTOSIS) == 33
        assert len(set(pathways.FERROPTOSIS)) == 33

    def test_includes_analysis_specific_genes(self) -> None:
        for gene in self.ANALYSIS_ONLY:
            assert gene in pathways.FERROPTOSIS, f"{gene} missing from ferroptosis set"

    def test_excludes_removed_genes(self) -> None:
        for gene in self.REMOVED:
            assert gene not in pathways.FERROPTOSIS, f"{gene} should have been removed"

    def test_core_regulators_retained(self) -> None:
        # Named in Methods 2.4; present in both sets, must survive the swap.
        for gene in ("GPX4", "ACSL4", "SLC7A11", "FTH1", "NFE2L2"):
            assert gene in pathways.FERROPTOSIS


class TestFerroptosisCanonicalValue:
    """Pin the chordoma ferroptosis CV the whole manuscript uses.

    The chordoma CSC ferroptosis CV was reported two ways: 5.61 (the analysis /
    Colab set, printed on Figure 2, the scVI decomposition, and gene-set
    sensitivity) and 5.387 (an earlier pan-cancer run that used a different
    ferroptosis gene set). Once the package is standardized to the analysis set,
    5.61 is canonical everywhere, and the pan-cancer table carries it so Figure 3
    and the text agree.

    Rank consequence: at 5.61 chordoma is seventh of twelve for ferroptosis (it was
    quoted as fifth under the old value), still moderate coordination, well inside
    the range of solid tumors.
    """

    def test_chordoma_ferroptosis_cv_is_canonical(self) -> None:
        if not DATA_CSV.exists():
            pytest.skip(f"{DATA_CSV} not present")
        panel = pd.read_csv(DATA_CSV).set_index("cancer_type")
        assert panel.loc["Chordoma", "ferroptosis_cv"] == pytest.approx(5.609, abs=1e-3)

    def test_chordoma_ferroptosis_ranks_seventh(self) -> None:
        # Results 3.4 (revised): "seventh of twelve (CV = 5.61)".
        if not DATA_CSV.exists():
            pytest.skip(f"{DATA_CSV} not present")
        panel = pd.read_csv(DATA_CSV).sort_values("ferroptosis_cv").reset_index(drop=True)
        rank = int(panel.index[panel["cancer_type"] == "Chordoma"][0]) + 1
        assert rank == 7

    def test_pan_cancer_ferroptosis_all_finite(self) -> None:
        # Every cancer type must carry a ferroptosis CV computed on the same
        # convention; a missing value would mean a row escaped the shared pipeline.
        if not DATA_CSV.exists():
            pytest.skip(f"{DATA_CSV} not present")
        panel = pd.read_csv(DATA_CSV)
        assert panel["ferroptosis_cv"].notna().all()
