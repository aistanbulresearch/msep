"""Check the assembled supplementary tables against the manuscript.

Every table the manuscript declares has a stated shape -- 42 pairwise comparisons,
23 defense genes, 12 cancer types, four gene sets of 33/29/28/20 genes. These tests
assert the assembled tables still have those shapes and still carry the values the
Results quote, so a rebuild cannot quietly disagree with the text.

The tests build the tables in-process rather than reading the emitted files, so they
run without ``figures/`` present. They skip when ``processed/`` is absent, since it
is not tracked in git.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "notebooks" / "paper_figures"))

if not (REPO_ROOT / "processed").exists():
    pytest.skip("processed/ not available (untracked)", allow_module_level=True)

import build_supplementary_tables as builder  # noqa: E402


@pytest.fixture(scope="module")
def tables() -> dict[str, pd.DataFrame]:
    return builder.build_all()


class TestDeclaredShapes:
    """Each table must have the size the manuscript says it has."""

    @pytest.mark.parametrize(
        "table,expected_rows",
        [
            ("S1", 110),   # 33 + 29 + 28 + 20 pathway genes
            ("S2", 10),    # "Ten cell types were identified in total" (3.1)
            ("S3", 10),    # per-cell entropy statistics for all ten cell types
            ("S4", 42),    # "all 42 pairwise Mann-Whitney U tests" (3.2)
            ("S5", 12),    # 12 cancer types (3.4)
            ("S6", 23),    # "the 23 key defense genes examined" (3.3)
            ("S7", 6),     # six key defense genes (2.7)
            ("S8", 9),     # nine cancer types (3.6)
            ("S10", 9),    # CSC vs three cell types across three pathways
        ],
    )
    def test_row_count(self, tables, table: str, expected_rows: int) -> None:
        assert len(tables[table]) == expected_rows

    def test_all_ten_tables_present(self, tables) -> None:
        assert sorted(tables, key=lambda s: int(s[1:])) == [
            f"S{i}" for i in range(1, 11)
        ]

    def test_no_table_is_empty(self, tables) -> None:
        for name, df in tables.items():
            assert not df.empty, f"Table {name} is empty"


class TestS1GeneSets:
    """Methods 2.4 states the size of each curated gene set."""

    @pytest.mark.parametrize(
        "pathway,expected",
        [("ferroptosis", 33), ("immune_evasion", 29), ("emt", 28), ("housekeeping", 20)],
    )
    def test_gene_set_size(self, tables, pathway: str, expected: int) -> None:
        assert (tables["S1"]["pathway"] == pathway).sum() == expected

    def test_genes_unique_within_pathway(self, tables) -> None:
        for _, group in tables["S1"].groupby("pathway"):
            assert group["gene"].is_unique

    def test_key_genes_present(self, tables) -> None:
        # Genes the Results discuss by name must actually be in the sets.
        s1 = tables["S1"]
        for gene, pathway in [("TBXT", "emt"), ("VIM", "emt"), ("GPX4", "ferroptosis"),
                              ("FTH1", "ferroptosis"), ("HLA-E", "immune_evasion"),
                              ("B2M", "immune_evasion")]:
            match = s1[(s1["gene"] == gene) & (s1["pathway"] == pathway)]
            assert len(match) == 1, f"{gene} missing from {pathway}"


class TestS2Annotation:
    def test_cell_counts_sum_to_cohort_total(self, tables) -> None:
        # "106,584 cells were retained" (3.1).
        assert tables["S2"]["n_cells"].sum() == 106_584

    def test_csc_is_the_analysis_population(self, tables) -> None:
        # "all entropy and CV analyses [...] were performed on the CSC_TBXT+
        # population (n=6,730)" (2.3).
        s2 = tables["S2"].set_index("cell_type")
        assert s2.loc["CSC_TBXT+", "n_cells"] == 6_730

    def test_malignant_compartment_total(self, tables) -> None:
        # "The broader malignant population [...] totaled 9,215 cells" (2.3).
        s2 = tables["S2"].set_index("cell_type")
        assert s2.loc["CSC_TBXT+", "n_cells"] + s2.loc["Malignant_TBXT-", "n_cells"] == 9_215


class TestS3EntropyStatistics:
    """Results 3.2 quotes per-cell entropy values for CSC_TBXT+."""

    @pytest.mark.parametrize(
        "column,quoted",
        [
            ("global_median", 9.97),
            ("no_ribo_median", 10.03),
            ("immune_evasion_norm_median", 0.409),
            ("ferroptosis_norm_median", 0.438),
            ("emt_mechanical_norm_median", 0.437),
            ("housekeeping_norm_median", 0.683),
        ],
    )
    def test_quoted_values_are_medians(self, tables, column: str, quoted: float) -> None:
        # Results 3.2 introduces 9.97 bits as "a mean". Every quoted value in that
        # section matches the median, not the mean -- the CSC global mean is 9.838.
        # Figure 2A already says "ranked by median entropy"; the prose needs the
        # same word.
        csc = tables["S3"].set_index("cell_type_fine").loc["CSC_TBXT+"]
        assert csc[column] == pytest.approx(quoted, abs=5e-3)

    def test_csc_has_highest_global_entropy(self, tables) -> None:
        # "CSC_TBXT+ cells showed the highest global entropy among all ten cell
        # types" (3.2).
        s3 = tables["S3"]
        assert s3.loc[s3["global_median"].idxmax(), "cell_type_fine"] == "CSC_TBXT+"

    def test_ribosomal_exclusion_preserves_the_ranking(self, tables) -> None:
        # "Ribosomal gene exclusion did not alter this ranking" (3.2).
        s3 = tables["S3"]
        assert s3.loc[s3["no_ribo_median"].idxmax(), "cell_type_fine"] == "CSC_TBXT+"

    def test_immune_evasion_entropy_is_lowest_pathway_in_csc(self, tables) -> None:
        # "CSC showed the lowest normalized per-cell entropy for immune evasion
        # genes (0.409)" (3.2).
        csc = tables["S3"].set_index("cell_type_fine").loc["CSC_TBXT+"]
        pathway_cols = [
            "ferroptosis_norm_median", "immune_evasion_norm_median",
            "emt_mechanical_norm_median", "housekeeping_norm_median",
        ]
        values = {c: csc[c] for c in pathway_cols}
        assert min(values, key=values.get) == "immune_evasion_norm_median"


class TestS4StatisticalTests:
    def test_all_comparisons_significant_after_fdr(self, tables) -> None:
        # "significantly higher than every other cell type in all 42 pairwise
        # Mann-Whitney U tests" (3.2).
        s4 = tables["S4"]
        assert (s4["p_adjusted"] < 0.05).all()


class TestS8Consolidation:
    def test_five_cancers_show_full_consolidation(self, tables) -> None:
        # "Full consolidation (3 of 3 pathways) was observed in five cancer types"
        # (3.6).
        s8 = tables["S8"]
        assert (s8["n_consolidated"] == 3).sum() == 5
