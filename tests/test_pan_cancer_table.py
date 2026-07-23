"""Tie the committed pan-cancer table to the claims made in the manuscript.

Figure 3 and Results 3.4 of Cavus & Kuskucu (2026) state specific CV values and
specific ranks for chordoma within the 12-cancer panel. These tests assert that
``data/pan_cancer_cv.csv`` still supports every one of those statements, so that
extending the panel (WP-3) or regenerating the table cannot silently drift away
from the text.

If a test here fails, either the table or the manuscript sentence it names has
to change -- they are no longer describing the same data.
"""

from pathlib import Path

import pandas as pd
import pytest

DATA_CSV = Path(__file__).resolve().parents[1] / "data" / "pan_cancer_cv.csv"

# Tolerance for values the manuscript quotes rounded to 2 decimals (e.g. "5.49").
ROUNDED_TOL = 0.005


@pytest.fixture(scope="module")
def panel() -> pd.DataFrame:
    """The committed pan-cancer table."""
    if not DATA_CSV.exists():
        pytest.skip(
            f"{DATA_CSV} not present; run notebooks/pan_cancer/build_pan_cancer_table.py"
        )
    return pd.read_csv(DATA_CSV)


def _rank_of(panel: pd.DataFrame, cancer_type: str, column: str) -> int:
    """1-based rank of ``cancer_type`` when ``column`` is sorted ascending."""
    ordered = panel.sort_values(column).reset_index(drop=True)
    return int(ordered.index[ordered["cancer_type"] == cancer_type][0]) + 1


class TestSchema:
    """The table must be structurally usable before its numbers mean anything."""

    def test_twelve_cancer_types(self, panel: pd.DataFrame) -> None:
        # Results 3.4: "malignant cells from 12 cancer types".
        assert len(panel) == 12

    def test_cancer_types_unique(self, panel: pd.DataFrame) -> None:
        assert panel["cancer_type"].is_unique

    def test_no_missing_cv_values(self, panel: pd.DataFrame) -> None:
        cv_cols = [c for c in panel.columns if c.endswith("_cv")]
        assert panel[cv_cols].notna().all().all()

    def test_every_row_has_provenance(self, panel: pd.DataFrame) -> None:
        for col in ("source", "accession", "platform", "cell_selection"):
            assert panel[col].notna().all(), f"{col} has gaps"

    def test_only_smartseq2_excluded_from_primary(self, panel: pd.DataFrame) -> None:
        # Methods 2.1.3: the primary comparison is 10x-only; Smart-seq2 melanoma
        # is reported separately because of platform effects on variance.
        excluded = panel.loc[~panel["in_primary_comparison"], "platform"]
        assert set(excluded) == {"Smart-seq2"}


class TestEmtCoordination:
    """Results 3.4, EMT paragraph."""

    def test_chordoma_lowest_emt_among_10x(self, panel: pd.DataFrame) -> None:
        # "Among the 10x-platform cancers, chordoma CSC exhibited the lowest
        # EMT CV (4.632)".
        primary = panel[panel["in_primary_comparison"]]
        winner = primary.loc[primary["emt_cv"].idxmin()]
        assert winner["cancer_type"] == "Chordoma"
        assert winner["emt_cv"] == pytest.approx(4.632, abs=1e-3)

    def test_next_three_lowest_emt_are_osteo_bcc_breast(self, panel: pd.DataFrame) -> None:
        # "...followed by basal cell carcinoma (5.49), osteosarcoma (5.25), and
        # breast cancer (5.49)". The manuscript lists these three out of
        # ascending order; the set and the values are what matter.
        primary = panel[panel["in_primary_comparison"]].sort_values("emt_cv")
        runners_up = primary.iloc[1:4]
        assert set(runners_up["cancer_type"]) == {"Osteosarcoma*", "BCC", "Breast"}
        expected = {"Osteosarcoma*": 5.25, "BCC": 5.49, "Breast": 5.49}
        for _, row in runners_up.iterrows():
            assert row["emt_cv"] == pytest.approx(
                expected[row["cancer_type"]], abs=ROUNDED_TOL
            )

    @pytest.mark.parametrize(
        "cancer_type,expected_cv",
        [("GBM", 6.72), ("Breast_TNBC", 7.64), ("Breast_Ductal", 7.83)],
    )
    def test_high_emt_heterogeneity_values(
        self, panel: pd.DataFrame, cancer_type: str, expected_cv: float
    ) -> None:
        # "Glioblastoma (6.72), triple-negative breast cancer (7.64), and
        # invasive ductal breast carcinoma (7.83) showed considerably higher
        # EMT heterogeneity".
        row = panel.loc[panel["cancer_type"] == cancer_type].iloc[0]
        assert row["emt_cv"] == pytest.approx(expected_cv, abs=ROUNDED_TOL)


class TestPathwaySelectivity:
    """Results 3.4: chordoma is exceptional for EMT but not for the other two."""

    def test_chordoma_ferroptosis_ranks_seventh_of_twelve(self, panel: pd.DataFrame) -> None:
        # "For ferroptosis, chordoma ranked seventh of twelve (CV = 5.61),
        # indicating moderate but not exceptional coordination." Value and rank
        # follow from standardizing on the analysis ferroptosis gene set.
        assert _rank_of(panel, "Chordoma", "ferroptosis_cv") == 7
        chordoma = panel.loc[panel["cancer_type"] == "Chordoma"].iloc[0]
        assert chordoma["ferroptosis_cv"] == pytest.approx(5.609, abs=1e-3)

    def test_chordoma_immune_evasion_ranks_ninth_of_twelve(
        self, panel: pd.DataFrame
    ) -> None:
        # "For immune evasion, chordoma ranked ninth of twelve (CV = 8.548),
        # placing it among the most heterogeneous cancer types".
        assert _rank_of(panel, "Chordoma", "immune_evasion_cv") == 9
        chordoma = panel.loc[panel["cancer_type"] == "Chordoma"].iloc[0]
        assert chordoma["immune_evasion_cv"] == pytest.approx(8.548, abs=1e-3)

    @pytest.mark.parametrize(
        "cancer_type,expected_cv",
        [("RCC_Nonpap", 8.57), ("Breast_TNBC", 9.33), ("GBM", 9.52)],
    )
    def test_immune_evasion_neighbours(
        self, panel: pd.DataFrame, cancer_type: str, expected_cv: float
    ) -> None:
        # "...alongside renal cell carcinoma (8.57), triple-negative breast
        # cancer (9.33), and glioblastoma (9.52)."
        row = panel.loc[panel["cancer_type"] == cancer_type].iloc[0]
        assert row["immune_evasion_cv"] == pytest.approx(expected_cv, abs=ROUNDED_TOL)

    def test_chordoma_is_selectively_coordinated(self, panel: pd.DataFrame) -> None:
        # The central pathway-selectivity claim: chordoma's EMT rank must be
        # better than both its ferroptosis and its immune evasion rank.
        emt_rank = _rank_of(panel, "Chordoma", "emt_cv")
        assert emt_rank < _rank_of(panel, "Chordoma", "ferroptosis_cv")
        assert emt_rank < _rank_of(panel, "Chordoma", "immune_evasion_cv")
