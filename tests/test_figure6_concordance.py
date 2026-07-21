"""Guard the statistics reported on manuscript Figure 6.

The version of panel D carried in the manuscript reported ``Spearman rho = 0.000,
p = 1.0000``. That figure was drawn from an 8-gene subset whose own correlation was
-0.500, while the full 23-gene source data gives -0.110. Nothing in the repository
would have caught the mismatch.

These tests recompute each statistic Figure 6 prints, straight from the analysis
output, so a printed value can no longer drift away from the data behind it.

The tests skip when ``processed/`` is absent, since it is not tracked in git.
"""

from pathlib import Path

import pandas as pd
import pytest
from scipy.stats import spearmanr

PROCESSED = Path(__file__).resolve().parents[1] / "processed"

# Values quoted in Results 3.7, rounded to 3 decimals in the manuscript.
QUOTED_TOL = 5e-4


def _read(name: str) -> pd.DataFrame:
    path = PROCESSED / name
    if not path.exists():
        pytest.skip(f"{path} not available (processed/ is untracked)")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def gene_concordance() -> pd.DataFrame:
    return _read("gene_cv_concordance_scRNA_vs_bulk.csv")


class TestSupplementaryGeneLevel:
    """Gene-level CV concordance -- Supplementary Figure S15.

    This analysis was formerly Figure 6 panel D. Recomputed over the full gene set
    the correlation is weak and not significant, so it no longer supports a
    main-figure concordance claim and now stands in the supplement as a reported
    null result. The pathway-level concordance in Figure 6A is unaffected.
    """

    def test_uses_all_twenty_three_genes(self, gene_concordance: pd.DataFrame) -> None:
        # Results 3.7 and Table S6 both describe 23 key defense genes. S15 must plot
        # the whole set, not a subset chosen after the fact.
        assert len(gene_concordance) == 23
        assert gene_concordance["gene"].is_unique

    def test_reported_correlation_matches_the_data(
        self, gene_concordance: pd.DataFrame
    ) -> None:
        rho, pval = spearmanr(
            gene_concordance["scRNA_CV"], gene_concordance["Bulk_CV"]
        )
        assert rho == pytest.approx(-0.110, abs=1e-3)
        assert pval == pytest.approx(0.618, abs=1e-3)

    def test_gene_level_concordance_is_not_significant(
        self, gene_concordance: pd.DataFrame
    ) -> None:
        # The honest result. Any narrative asserting gene-level cross-platform
        # concordance has to be reconciled with this, not the other way round.
        _, pval = spearmanr(
            gene_concordance["scRNA_CV"], gene_concordance["Bulk_CV"]
        )
        assert pval >= 0.05

    def test_slc7a11_is_the_most_variable_gene_in_scrna(
        self, gene_concordance: pd.DataFrame
    ) -> None:
        # Results 3.7 groups SLC7A11 with genes whose bulk CV is low "consistent with
        # their low CV in the scRNA-seq data". Its scRNA-seq CV is in fact the highest
        # of all 23. This test pins the fact so the sentence gets revisited.
        top = gene_concordance.loc[gene_concordance["scRNA_CV"].idxmax()]
        assert top["gene"] == "SLC7A11"
        assert top["scRNA_CV"] == pytest.approx(13.73, abs=1e-2)


class TestPanelA:
    """Cross-platform pathway CV concordance."""

    @pytest.mark.parametrize(
        "pathway,expected",
        [
            ("emt_mechanical", 0.784),
            ("housekeeping", 0.891),
            ("ferroptosis", 0.904),
            ("immune_evasion", 1.113),
        ],
    )
    def test_bulk_pathway_cv_matches_manuscript(
        self, pathway: str, expected: float
    ) -> None:
        # "EMT showed the lowest across-tumor CV (0.784), followed by housekeeping
        # (0.891), ferroptosis (0.904), and immune evasion (1.113)."
        bulk = _read("bulk_validation_cv.csv").set_index("group")
        assert bulk.loc["Chordoma", f"{pathway}_CV"] == pytest.approx(
            expected, abs=QUOTED_TOL
        )

    def test_pathway_hierarchy_is_preserved_across_platforms(self) -> None:
        # The load-bearing claim of section 3.7: EMT most homogeneous, immune evasion
        # most heterogeneous, in bulk exactly as in scRNA-seq.
        bulk = _read("bulk_validation_cv.csv").set_index("group")
        resistance = ["emt_mechanical", "ferroptosis", "immune_evasion"]
        ordered = sorted(resistance, key=lambda p: bulk.loc["Chordoma", f"{p}_CV"])
        assert ordered == ["emt_mechanical", "ferroptosis", "immune_evasion"]

        pan = _read("pan_cancer_census_cv_expanded.csv").set_index("cancer_type")
        ordered_sc = sorted(resistance, key=lambda p: pan.loc["Chordoma", f"{p}_CV"])
        assert ordered_sc == ordered


class TestPanelB:
    """Chordoma vs notochord."""

    def test_notochord_emt_cv_is_lowest_among_resistance_pathways(self) -> None:
        # Results 3.7 says "notochord EMT CV (0.282) was the lowest among all
        # notochord pathways". Taken literally that is wrong: housekeeping is lower
        # (0.207). EMT is the lowest of the three *resistance* pathways, which is the
        # comparison the surrounding argument actually needs -- the manuscript uses
        # the same resistance-only framing for chordoma in section 3.3, where
        # housekeeping is described as an included reference rather than a pathway
        # under test. The sentence needs the qualifier; the biology is unaffected.
        bulk = _read("bulk_validation_cv.csv").set_index("group")
        resistance = ["emt_mechanical", "ferroptosis", "immune_evasion"]
        values = {p: bulk.loc["Notochord", f"{p}_CV"] for p in resistance}
        assert min(values, key=values.get) == "emt_mechanical"
        assert values["emt_mechanical"] == pytest.approx(0.282, abs=QUOTED_TOL)

    def test_notochord_housekeeping_is_lower_than_emt(self) -> None:
        # The fact that makes the unqualified claim above inaccurate.
        bulk = _read("bulk_validation_cv.csv").set_index("group")
        assert (
            bulk.loc["Notochord", "housekeeping_CV"]
            < bulk.loc["Notochord", "emt_mechanical_CV"]
        )


class TestPanelC:
    """Per-patient stability."""

    def test_emt_most_homogeneous_in_four_of_six_patients(self) -> None:
        # "EMT was the most homogeneous resistance pathway in four of six patients
        # (NU02854, NU03174, NU03231, NU03372)."
        per_patient = _read("per_patient_pathway_cv.csv")
        emt_first = set(
            per_patient.loc[
                per_patient["most_homogeneous"] == "emt_mechanical", "patient_id"
            ]
        )
        assert emt_first == {"NU02854", "NU03174", "NU03231", "NU03372"}

    def test_nu02757_ferroptosis_marginally_beats_emt(self) -> None:
        # "In patient NU02757, ferroptosis was marginally more homogeneous than EMT
        # (CV 5.22 versus 5.30, a difference of 0.08)."
        per_patient = _read("per_patient_pathway_cv.csv").set_index("patient_id")
        row = per_patient.loc["NU02757"]
        assert row["ferroptosis_CV"] == pytest.approx(5.22, abs=0.005)
        assert row["emt_mechanical_CV"] == pytest.approx(5.30, abs=0.005)
        assert row["emt_mechanical_CV"] - row["ferroptosis_CV"] == pytest.approx(
            0.08, abs=0.005
        )
