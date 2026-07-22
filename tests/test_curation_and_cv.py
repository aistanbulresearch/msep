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


class TestFerroptosisCanonicalValue:
    """Pin the one ferroptosis CV the manuscript should use everywhere.

    Chordoma CSC ferroptosis CV appears in the analysis output as two values:
    5.387 (used by the text section 3.3/3.4 and by the pan-cancer table) and 5.609
    (in the scVI and gene-set-sensitivity intermediates, and printed on the
    published Figure 2). The two come from runs that resolved a different number of
    ferroptosis genes for the chordoma cells -- one found the whole 33-symbol set,
    the other 32.

    The pan-cancer comparison and the two rank claims in section 3.4 all rest on
    5.387 computed the same way across every cancer type, so 5.387 is canonical.
    Figure 2's 5.61 is the outlier and must be brought to 5.39. This test fixes the
    canonical value in the tracked data so the rest cannot drift away from it.
    """

    def test_chordoma_ferroptosis_cv_is_canonical(self) -> None:
        if not DATA_CSV.exists():
            pytest.skip(f"{DATA_CSV} not present")
        panel = pd.read_csv(DATA_CSV).set_index("cancer_type")
        assert panel.loc["Chordoma", "ferroptosis_cv"] == pytest.approx(5.387, abs=1e-3)

    def test_pan_cancer_ferroptosis_all_finite(self) -> None:
        # Every cancer type must carry a ferroptosis CV computed on the same
        # convention; a missing value would mean a row escaped the shared pipeline.
        if not DATA_CSV.exists():
            pytest.skip(f"{DATA_CSV} not present")
        panel = pd.read_csv(DATA_CSV)
        assert panel["ferroptosis_cv"].notna().all()
