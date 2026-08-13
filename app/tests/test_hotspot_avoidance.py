import pandas as pd
import pytest

from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult
from modules.hotspot_avoidance.hotspot_avoidance_main import make_dcub_custom_score
from modules.hotspot_avoidance.hotspot_avoidance_main import patch_sequence
from modules.hotspot_avoidance.hotspot_avoidance_main import widen_slippage_base_units
from modules.shared_functions_and_vars import synonymous_codons
from modules.shared_functions_and_vars import translate

#: The exact column set eso.detection.dispatch.find_slippage_sites emits.
SLIPPAGE_COLUMNS = ["start", "end", "length_base_unit", "sequence", "num_base_units",
                    "log10_prob_slippage_ecoli"]


def _slippage_df(*rows):
    return pd.DataFrame(list(rows), columns=SLIPPAGE_COLUMNS)


# The flanks of the planted-hotspot sequence. Both are 30nt / 10 codons and,
# crucially, are themselves free of every site ESO's detectors look for - so
# the ONLY hotspot window in `prefix + "AAA"*5 + suffix` is the planted
# homopolymer at [30, 45). Without that, ESO's PSSM-based motif detector fires
# on near-matches (a bare "GATG" scores positive against the dam GATC motif),
# those flank positions become legitimately editable, and the locality
# assertion below would be testing nothing.
PLANTED_PREFIX = "ATGCCACAACACGCACGCAGCTACAACGTG"
PLANTED_SUFFIX = "CAAGTCTCACTAGTGAGTGACTTCGGTAAT"

# 20 distinct codons with no repeats, no homopolymer runs, and - verified
# against the detectors themselves - no recombination, slippage, or motif
# hits at all, so `patch_sequence` takes the "nothing detected" early return.
CLEAN_SEQUENCE = "ATGGTTACTAGTTGTCCGGCAACCAACGTAACGAGAGGCTATGTGCAAGCCGGAATTCTC"


def _flat_codon_table(preferred=()):
    """A table where every codon scores 1.0 except a few preferred ones at 2.0 -
    enough to give the optimizer a gradient without needing real organisms."""
    table = {
        amino_acid: {codon: 1.0 for codon in codons}
        for amino_acid, codons in synonymous_codons.items()
    }
    for codon in preferred:
        from modules.shared_functions_and_vars import nt_to_aa

        table[nt_to_aa[codon]][codon] = 2.0
    return table


class TestMakeDcubCustomScore:
    def test_sums_per_codon_table_lookups(self):
        score = make_dcub_custom_score(_flat_codon_table(preferred=["TGT"]))
        # TGT (2.0) + AAA (1.0) + ATG (1.0)
        assert score("TGTAAAATG") == pytest.approx(4.0)

    def test_unknown_codon_contributes_zero(self):
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("NNNAAA") == pytest.approx(1.0)

    def test_trailing_partial_codon_is_ignored(self):
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("AAAAA") == pytest.approx(1.0)

    def test_empty_sequence_scores_zero(self):
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("") == pytest.approx(0.0)


class TestWidenSlippageBaseUnits:
    """ESO's `exclusion_site_correcter` drops every avoidance row narrower than
    3nt as soon as any exclusion region is passed, and this module always passes
    them. Widening sub-codon base units is what keeps homopolymer repair alive."""

    def test_homopolymer_base_unit_widens_from_one_to_three(self):
        sequence = PLANTED_PREFIX + "AAA" * 5 + PLANTED_SUFFIX
        # The `sequence` field is deliberately wrong here: if the widened row
        # carried it over instead of recomputing it from the coordinates, the
        # marker would survive and the last assertion would catch it.
        df, warnings_out = widen_slippage_base_units(
            _slippage_df([30, 45, 1, "X" * 15, 15, -1.965]),
            sequence,
        )

        assert warnings_out == []
        assert len(df) == 1
        row = df.iloc[0]
        assert int(row["length_base_unit"]) == 3
        assert int(row["num_base_units"]) == 5
        assert int(row["start"]) == 30
        assert int(row["end"]) == 45
        assert row["sequence"] == "AAAAAAAAAAAAAAA" == sequence[30:45]
        # The column survives (optimization_engine indexes by name) but is
        # blanked: ESO derives it from the base unit and repeat count, so the
        # detected value would no longer satisfy that formula on a
        # re-parameterized row.
        assert list(df.columns) == SLIPPAGE_COLUMNS
        assert pd.isna(row["log10_prob_slippage_ecoli"])

    def test_dinucleotide_base_unit_widens_from_two_to_six(self):
        # These coordinates are real find_slippage_sites output for this sequence.
        sequence = PLANTED_PREFIX + "TG" * 7 + "CA" + PLANTED_SUFFIX
        df, warnings_out = widen_slippage_base_units(
            _slippage_df([27, 43, 2, "GTGTGTGTGTGTGTGT", 8, -4.245]),
            sequence,
        )

        assert warnings_out == []
        row = df.iloc[0]
        assert int(row["length_base_unit"]) == 6
        # (43 - 27) // 6 == 2, so the row shortens to cover 2 whole units.
        assert int(row["num_base_units"]) == 2
        assert int(row["start"]) == 27
        assert int(row["end"]) == 39
        assert row["sequence"] == sequence[27:39]
        assert len(row["sequence"]) == 12

    def test_base_unit_of_four_passes_through_untouched(self):
        """Regression: widening anything already >= 3nt breaks repairs ESO
        would have made. A 4nt unit clears ESO's `start < end - 2` filter as
        detected, but widening it to lcm(4, 3) == 12 leaves (42 - 30) // 12 == 1
        unit - below the 2-unit minimum - so the site would be dropped with a
        warning instead of repaired."""
        sequence = PLANTED_PREFIX + "ACGT" * 3 + PLANTED_SUFFIX
        original = _slippage_df([30, 42, 4, "ACGTACGTACGT", 3, -4.56])
        df, warnings_out = widen_slippage_base_units(original, sequence)

        assert warnings_out == []
        assert df.iloc[0].to_dict() == original.iloc[0].to_dict()
        assert int(df.iloc[0]["length_base_unit"]) == 4
        assert df.iloc[0]["log10_prob_slippage_ecoli"] == pytest.approx(-4.56)

    def test_codon_width_base_unit_passes_through_untouched(self):
        sequence = PLANTED_PREFIX + "GCT" * 5 + PLANTED_SUFFIX
        # A marker rather than the true slice: passing through leaves it alone,
        # whereas the widening branch would overwrite it with sequence[30:45].
        # Without this the test would pass either way, since lcm(3, 3) == 3.
        original = _slippage_df([30, 45, 3, "MARKER-NOT-A-SLICE", 5, -2.5])
        df, warnings_out = widen_slippage_base_units(original, sequence)

        assert warnings_out == []
        assert df.iloc[0].to_dict() == original.iloc[0].to_dict()
        assert df.iloc[0]["sequence"] == "MARKER-NOT-A-SLICE"

    def test_run_too_short_to_widen_is_dropped_with_a_warning(self):
        sequence = PLANTED_PREFIX + "AAAA" + PLANTED_SUFFIX
        # 4nt at 1nt base units -> (34 - 30) // 3 == 1 widened unit, and
        # modify_df_slippage needs at least 2 to emit anything.
        df, warnings_out = widen_slippage_base_units(
            _slippage_df([30, 34, 1, "AAAA", 4, -0.5]),
            sequence,
        )

        assert df.empty
        assert list(df.columns) == SLIPPAGE_COLUMNS
        assert len(warnings_out) == 1
        assert "30-34" in warnings_out[0]
        assert "1nt" in warnings_out[0]
        assert "left unmodified" in warnings_out[0]

    def test_empty_and_missing_dataframes_are_passed_straight_back(self):
        empty = _slippage_df()
        df, warnings_out = widen_slippage_base_units(empty, "ATGATG")
        assert df is empty
        assert warnings_out == []

        df, warnings_out = widen_slippage_base_units(None, "ATGATG")
        assert df is None
        assert warnings_out == []


class TestPatchSequence:
    def test_planted_homopolymer_is_removed_and_everything_else_is_untouched(self):
        """The core locality guarantee: the planted slippage hotspot is gone,
        and every nucleotide outside its widened window is byte-identical."""
        # 10 varied codons, then a 5x AAA homopolymer run (a slippage hotspot),
        # then 10 more varied codons.
        prefix = PLANTED_PREFIX                          # 10 codons, 30nt
        hotspot = "AAA" * 5                              # 15nt homopolymer run
        suffix = PLANTED_SUFFIX                          # 10 codons, 30nt
        sequence = prefix + hotspot + suffix
        assert len(sequence) % 3 == 0

        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        assert isinstance(result, HotspotPatchResult)
        assert result.sequence_before == sequence
        assert len(result.sequence_after) == len(sequence)
        # Translation is preserved - this is what makes an index-wise diff valid.
        assert translate(result.sequence_after) == translate(sequence)
        # The homopolymer run is broken.
        assert result.num_edits > 0
        assert "AAAAAAAAAAAAAAA" not in result.sequence_after
        # Everything outside the hotspot window is locked. The window is the
        # homopolymer widened to codon boundaries, i.e. [30, 45).
        assert result.sequence_after[:30] == prefix
        assert result.sequence_after[45:] == suffix

    def test_clean_sequence_is_returned_unchanged_with_zero_edits(self):
        sequence = CLEAN_SEQUENCE
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert result.sequence_after == sequence
        assert result.num_edits == 0

    def test_initiation_prefix_is_never_modified(self):
        """A hotspot overlapping the initiation-optimized prefix must not make
        that prefix editable - InitiationModule optimized it for weak folding."""
        prefix = "AAA" * 15                              # 45nt, all locked, and a hotspot
        suffix = "ATGGCTTGTGATGAACATATCAAGCTGAAT"
        sequence = prefix + suffix

        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=15,
        )
        assert result.sequence_after[:45] == prefix

    def test_result_reports_detected_site_counts(self):
        sequence = PLANTED_PREFIX + "AAA" * 5 + PLANTED_SUFFIX
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert set(result.detected_sites) == {"recombination", "slippage", "motifs"}
        assert result.detected_sites["slippage"] >= 1

    def test_summary_shape_matches_the_frontend_contract(self):
        sequence = "ATGGCTTGTGATGAACATATCAAGCTGAAT"
        summary = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        ).summary
        assert summary["enabled"] is True
        assert set(summary) == {
            "enabled",
            "sequence_before",
            "sequence_after",
            "num_edits",
            "detected_sites",
            "warnings",
        }

    def test_unwidenable_run_is_reported_but_still_counted_as_detected(self):
        """The drop path, end to end through patch_sequence rather than through
        the helper in isolation. A 2nt base unit widens to 6nt, and this run is
        only 8nt long - one whole unit, below modify_df_slippage's minimum of
        two - so it cannot be disrupted at codon resolution.

        Pins all three wiring behaviours at once: the warning reaches the user,
        `detected_sites` still reports the site as DETECTED (it is counted from
        the original detection, not from the widened rows), and the sequence
        comes back untouched rather than half-edited.
        """
        sequence = PLANTED_PREFIX + "CACACA" + PLANTED_SUFFIX
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        # Detected - and still reported as such even though it was not repaired.
        assert result.detected_sites["slippage"] == 1
        assert result.warnings == [
            "Slippage site at 30-38 (base unit 2nt) is too short to disrupt at "
            "codon resolution and was left unmodified."
        ]
        assert result.sequence_after == sequence
        assert result.num_edits == 0

    def test_four_nucleotide_base_unit_repeat_is_repaired(self):
        """Regression for the over-broad widening predicate: a 4nt base unit is
        already wide enough for ESO's exclusion filter, so it must be passed
        through as detected and actually repaired. Widening it to lcm(4, 3) == 12
        left only one whole unit, which silently downgraded a repairable site to
        an unrepairable one - with a warning that was itself wrong."""
        repeat = "ACGT" * 3
        sequence = PLANTED_PREFIX + repeat + PLANTED_SUFFIX
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        assert result.detected_sites["slippage"] == 1
        assert result.warnings == []
        assert result.num_edits > 0
        assert repeat not in result.sequence_after
        assert translate(result.sequence_after) == translate(sequence)

    def test_customscore_performance_warning_is_not_surfaced_to_users(self):
        """ESO's CustomScore always warns about per-trial rescoring cost. That is
        an implementation detail, not something a biologist should see next to
        genuine 'could not clear this hotspot' warnings."""
        sequence = "ATGGCTTGTGATGAACATATCAAGCTGAAT"
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert not any("re-evaluates score_fn" in warning for warning in result.warnings)
