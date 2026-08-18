import math

import pandas as pd
import pytest

from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult
from modules.hotspot_avoidance.hotspot_avoidance_main import make_dcub_custom_score
from modules.hotspot_avoidance.hotspot_avoidance_main import patch_sequence
from modules.hotspot_avoidance.hotspot_avoidance_main import widen_slippage_base_units
from modules.shared_functions_and_vars import nt_to_aa
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
    """Scores with general_geomean - DCUB's canonical routine - not by summing
    the table. Weights are the CAI-style (0, 1] scale loss_table_to_weights
    produces, on which a fully DCUB-optimal ORF scores 1.0."""

    def test_scores_the_geometric_mean_of_codon_weights(self):
        score = make_dcub_custom_score(_flat_codon_table(preferred=["TGT"]))
        # TGT (2.0) and AAA (1.0); ATG is non-synonymous, so general_geomean
        # skips it - there is no codon choice at a Met to optimize.
        assert score("TGTAAAATG") == pytest.approx(math.sqrt(2.0))

    def test_a_single_bad_codon_costs_more_than_it_would_by_summing(self):
        """The reason for the geometric mean: it is multiplicative, so one very
        poor codon drags the whole window down instead of being averaged away by
        its neighbours. This is also how EvaluationModule scores."""
        table = _flat_codon_table()
        table[nt_to_aa["TGT"]]["TGT"] = 0.01

        score = make_dcub_custom_score(table)
        all_good = score("AAAGATTTT")
        one_bad = score("TGTGATTTT")

        assert one_bad < all_good
        # gmean(0.01, 1, 1) = 0.01^(1/3) ~ 0.215, i.e. a ~78% drop from one
        # codon out of three. A sum would have dropped only 33%.
        assert one_bad == pytest.approx((0.01 * 1.0 * 1.0) ** (1 / 3))

    def test_unknown_codon_falls_back_to_the_mean_weight(self):
        """general_geomean substitutes the profile's mean rather than dropping
        the codon or scoring it zero - a zero would zero the whole sequence."""
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("NNNAAA") == pytest.approx(1.0)

    def test_trailing_partial_codon_is_ignored(self):
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("AAAAA") == pytest.approx(1.0)

    def test_too_short_to_hold_a_codon_scores_zero(self):
        """gmean of an empty list is nan, which would poison DNAChisel's
        objective comparison. ESO can hand in a slice shorter than a codon."""
        score = make_dcub_custom_score(_flat_codon_table())
        assert score("") == pytest.approx(0.0)
        assert score("AT") == pytest.approx(0.0)


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

    def test_run_shorter_than_one_widened_unit_is_dropped_with_a_warning(self):
        """The remaining unrepairable case: not even ONE widened unit fits, so
        there is no chunk that both covers whole codons and stays inside the
        detected repeat. 4nt at 2nt base units -> (34 - 30) // 6 == 0."""
        sequence = PLANTED_PREFIX + "ACAC" + PLANTED_SUFFIX
        df, warnings_out = widen_slippage_base_units(
            _slippage_df([30, 34, 2, "ACAC", 2, -0.5]),
            sequence,
        )

        assert df.empty
        assert list(df.columns) == SLIPPAGE_COLUMNS
        assert len(warnings_out) == 1
        assert "30-34" in warnings_out[0]
        assert "2nt" in warnings_out[0]
        assert "left unmodified" in warnings_out[0]

    def test_a_single_widened_unit_is_enough(self):
        """A run holding exactly ONE widened unit used to be discarded, because
        modify_df_slippage iterates range(0, num_base_units - 1, 2) and emits
        nothing below 2. Declaring 2 units makes that loop emit chunk 0 - which
        IS the real repeat - so the minimum repairable run halves: 6nt for a
        dinucleotide instead of 12nt."""
        sequence = PLANTED_PREFIX + "CGCGCG" + PLANTED_SUFFIX
        df, warnings_out = widen_slippage_base_units(
            _slippage_df([30, 36, 2, "CGCGCG", 3, -4.2]),
            sequence,
        )

        assert warnings_out == []
        assert len(df) == 1
        row = df.iloc[0]
        assert int(row["length_base_unit"]) == 6
        assert int(row["num_base_units"]) == 2
        # `end` stays at the real extent of the repeat - only the declared unit
        # COUNT is inflated, never the coordinates.
        assert int(row["start"]) == 30
        assert int(row["end"]) == 36
        assert row["sequence"] == "CGCGCG" == sequence[30:36]

    def test_declared_units_never_extend_the_avoided_chunk_past_the_repeat(self):
        """The property that makes declaring 2 units sound and padding the
        window unsound: every chunk ESO derives must lie INSIDE the detected
        repeat. A padded row would let DNAChisel satisfy the constraint by
        editing the flank and leaving the repeat intact - measured on mCherry,
        padding one codon each side around CGCGCG produced GAG->GAA, one
        reported edit, and the repeat fully intact."""
        from eso.detection.slippage import modify_df_slippage

        sequence = PLANTED_PREFIX + "CGCGCG" + PLANTED_SUFFIX
        df, _ = widen_slippage_base_units(
            _slippage_df([30, 36, 2, "CGCGCG", 3, -4.2]),
            sequence,
        )

        chunks = modify_df_slippage(df)
        assert len(chunks) == 1, "chunk 1 must never be emitted"
        chunk = chunks.iloc[0]
        assert (int(chunk["start"]), int(chunk["end"])) == (30, 36)
        assert chunk["sequence"] == "CGCGCG"
        # Wider than ESO's `df[df.start < df.end - 2]` filter, which is the
        # whole reason widening exists at all.
        assert int(chunk["start"]) < int(chunk["end"]) - 2

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
            "detected_regions",
            "warnings",
        }


class TestDetectedRegions:
    """The regions drive the frontend's highlight overlay, so their coordinate
    convention is load-bearing: 0-indexed, exclusive end, into sequence_before."""

    def test_regions_are_reported_with_kind_and_coordinates(self):
        # A 15x"A" homopolymer, in frame, with clean flanks.
        sequence = "ATG" + "AAA" * 5 + "GCTTGTGATGAACAT"
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        assert result.detected_regions, "the planted homopolymer must be reported"
        for region in result.detected_regions:
            assert set(region) == {"kind", "start", "end"}
            assert region["kind"] in {"recombination", "slippage", "motifs"}
            assert 0 <= region["start"] < region["end"] <= len(sequence)

    def test_region_slices_sequence_before_at_the_hotspot(self):
        """The decisive property: slicing sequence_before with a reported
        region must return the hotspot itself. An off-by-one or an inclusive
        end would put the highlight on the wrong nucleotides."""
        sequence = "ATG" + "AAA" * 5 + "GCTTGTGATGAACAT"
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        slippage = [r for r in result.detected_regions if r["kind"] == "slippage"]
        assert slippage, "a 15nt homopolymer is a slippage site"
        for region in slippage:
            assert set(result.sequence_before[region["start"]:region["end"]]) == {"A"}

    def test_counts_and_regions_agree(self):
        sequence = "ATG" + "AAA" * 5 + "GCTTGTGATGAACAT"
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        # Recombination reports a PAIR of windows per detected row, so its
        # region count is twice its site count; the other two are 1:1.
        kinds = [region["kind"] for region in result.detected_regions]
        assert kinds.count("slippage") == result.detected_sites["slippage"]
        assert kinds.count("motifs") == result.detected_sites["motifs"]
        assert kinds.count("recombination") == 2 * result.detected_sites["recombination"]

    def test_a_clean_sequence_reports_no_regions(self):
        result = patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert result.detected_regions == []
        assert result.summary["detected_regions"] == []

    def test_unwidenable_run_is_reported_but_still_counted_as_detected(self, monkeypatch):
        """The drop path, end to end through patch_sequence.

        Detection is stubbed rather than planted in the sequence because every
        sub-codon site ESO's detector actually emits is now repairable: its
        smallest reported span is 12nt for a 1nt base unit and 6nt for a 2nt
        one, both clearing the one-widened-unit minimum. The code path still
        exists for a site narrower than lcm(base_unit, 3), so it stays covered
        - it just can no longer be reached with real detector output.

        Pins all three wiring behaviours at once: the warning reaches the user,
        `detected_sites` still reports the site as DETECTED (it is counted from
        the original detection, not from the widened rows), and the sequence
        comes back untouched rather than half-edited.
        """
        from modules.hotspot_avoidance import hotspot_avoidance_main

        # 6nt filler keeps the ORF a whole number of codons; the stubbed row
        # below describes a 4nt window inside it.
        sequence = PLANTED_PREFIX + "ACACAC" + PLANTED_SUFFIX

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            return {
                "df_recombination": pd.DataFrame(),
                # 4nt at a 2nt base unit: (34 - 30) // lcm(2, 3) == 0, so not
                # even one widened unit fits.
                "df_slippage": _slippage_df([30, 34, 2, "ACAC", 2, -0.5]),
            }

        monkeypatch.setattr(
            hotspot_avoidance_main, "suspect_site_extractor", fake_extractor
        )

        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        # Detected - and still reported as such even though it was not repaired.
        assert result.detected_sites["slippage"] == 1
        assert result.warnings == [
            "Slippage site at 30-34 (base unit 2nt) is too short to disrupt at "
            "codon resolution and was left unmodified."
        ]
        assert result.sequence_after == sequence
        assert result.num_edits == 0

    def test_six_nucleotide_dinucleotide_repeat_is_now_repaired(self):
        """The improvement, end to end: a 6nt CG repeat used to need 12nt to be
        touchable at all. It is now repaired, and only inside its own window."""
        sequence = PLANTED_PREFIX + "CGCGCG" + PLANTED_SUFFIX
        result = patch_sequence(
            sequence=sequence,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        assert result.detected_sites["slippage"] >= 1
        assert result.num_edits > 0
        assert result.warnings == []
        assert translate(result.sequence_before) == translate(result.sequence_after)
        # Everything outside the codon-widened detected window is untouched.
        editable_start = (min(r["start"] for r in result.detected_regions) // 3) * 3
        editable_end = -(-max(r["end"] for r in result.detected_regions) // 3) * 3
        assert result.sequence_after[:editable_start] == sequence[:editable_start]
        assert result.sequence_after[editable_end:] == sequence[editable_end:]

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


class TestDetectionConfiguration:
    def test_motifs_are_off_by_default(self, monkeypatch):
        """Motif detection is off unless explicitly enabled - ESO's PSSM scan
        reports far more sites than are biologically real, and each one is a
        licence to rewrite DCUB's chosen codons."""
        captured = {}

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            captured["compute_motifs"] = compute_motifs
            captured["common_motifs"] = kwargs.get("common_motifs")
            import pandas as pd
            return {"df_recombination": pd.DataFrame(), "df_slippage": pd.DataFrame()}

        monkeypatch.setattr(
            "modules.hotspot_avoidance.hotspot_avoidance_main.suspect_site_extractor",
            fake_extractor,
        )
        patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert captured["compute_motifs"] is False

    def test_motifs_can_be_enabled_explicitly(self, monkeypatch):
        captured = {}

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            captured["compute_motifs"] = compute_motifs
            captured["common_motifs"] = kwargs.get("common_motifs")
            import pandas as pd
            return {
                "df_recombination": pd.DataFrame(),
                "df_slippage": pd.DataFrame(),
                "df_motifs": pd.DataFrame(),
            }

        monkeypatch.setattr(
            "modules.hotspot_avoidance.hotspot_avoidance_main.suspect_site_extractor",
            fake_extractor,
        )
        patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
            compute_motifs=True,
            common_motifs=["dam"],
        )
        assert captured["compute_motifs"] is True
        assert captured["common_motifs"] == ["dam"]

    def test_config_defaults_are_used_when_no_override_is_given(self, monkeypatch):
        """The config is the single source of truth for defaults, so an operator
        can enable motifs deployment-wide without a code change."""
        from modules.hotspot_avoidance import hotspot_avoidance_main

        captured = {}

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            captured["compute_motifs"] = compute_motifs
            captured["common_motifs"] = kwargs.get("common_motifs")
            import pandas as pd
            return {
                "df_recombination": pd.DataFrame(),
                "df_slippage": pd.DataFrame(),
                "df_motifs": pd.DataFrame(),
            }

        monkeypatch.setattr(
            hotspot_avoidance_main, "suspect_site_extractor", fake_extractor
        )
        monkeypatch.setitem(
            hotspot_avoidance_main.config["HOTSPOT_AVOIDANCE"], "COMPUTE_MOTIFS", True
        )
        monkeypatch.setitem(
            hotspot_avoidance_main.config["HOTSPOT_AVOIDANCE"], "COMMON_MOTIFS", ["dcm"]
        )
        patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert captured["compute_motifs"] is True
        assert captured["common_motifs"] == ["dcm"]

    def test_detection_modes_are_passed_through(self, monkeypatch):
        from modules.hotspot_avoidance import hotspot_avoidance_main

        captured = {}

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            captured.update(kwargs)
            import pandas as pd
            return {"df_recombination": pd.DataFrame(), "df_slippage": pd.DataFrame()}

        monkeypatch.setattr(
            hotspot_avoidance_main, "suspect_site_extractor", fake_extractor
        )
        patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
            recombination_mode="fast",
            slippage_mode="fast",
        )
        assert captured["recombination_mode"] == "fast"
        assert captured["slippage_mode"] == "fast"

    def test_detection_modes_are_omitted_when_not_overridden(self, monkeypatch):
        """Not passed at all, rather than passed as a value restated here - so
        ESO's own defaults apply. Restating them in configuration.yaml only
        created a second place to drift from if ESO ever changed them."""
        from modules.hotspot_avoidance import hotspot_avoidance_main

        captured = {}

        def fake_extractor(target_seq, compute_motifs, num_sites, **kwargs):
            captured.update(kwargs)
            captured["_keys"] = set(kwargs)
            import pandas as pd
            return {"df_recombination": pd.DataFrame(), "df_slippage": pd.DataFrame()}

        monkeypatch.setattr(
            hotspot_avoidance_main, "suspect_site_extractor", fake_extractor
        )
        patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )

        assert "recombination_mode" not in captured["_keys"]
        assert "slippage_mode" not in captured["_keys"]

    def test_enabling_motifs_without_configuring_any_is_rejected(self, monkeypatch):
        """COMMON_MOTIFS is NOT redundant with ESO's default. ESO guards with
        `if common_motifs:`, so its None default leaves the motif list empty and
        find_motif_sites returns an empty frame - motif detection would be
        switched on and silently do nothing. Fail loudly instead."""
        from modules.hotspot_avoidance import hotspot_avoidance_main

        monkeypatch.setitem(
            hotspot_avoidance_main.config["HOTSPOT_AVOIDANCE"], "COMPUTE_MOTIFS", True
        )
        monkeypatch.setitem(
            hotspot_avoidance_main.config["HOTSPOT_AVOIDANCE"], "COMMON_MOTIFS", []
        )

        with pytest.raises(ValueError, match="no motifs are configured"):
            patch_sequence(
                sequence=CLEAN_SEQUENCE,
                codon_table=_flat_codon_table(),
                skipped_codons_num=0,
            )

    def test_eso_default_common_motifs_would_detect_nothing(self):
        """Pins the upstream behaviour the guard above exists for, so that if a
        future ESO starts defaulting to its full bundled motif set, this test
        fails and the guard can be revisited rather than quietly outliving its
        reason."""
        import numpy as np
        from eso import suspect_site_extractor as real_extractor

        sequence = "ATG" + "GATC" * 5 + "CCAGG" * 4 + "GCTTGTGATGAACATATCAAG"

        without = real_extractor(sequence, compute_motifs=True, num_sites=np.inf)
        with_motifs = real_extractor(
            sequence, compute_motifs=True, num_sites=np.inf, common_motifs=["dam", "dcm"]
        )

        assert len(without["df_motifs"]) == 0
        assert len(with_motifs["df_motifs"]) > 0

    def test_detected_sites_reports_zero_motifs_when_disabled(self):
        """With motifs off, suspect_site_extractor omits df_motifs entirely -
        the count must still be present and zero, not missing, since the
        frontend contract requires all three keys."""
        result = patch_sequence(
            sequence=CLEAN_SEQUENCE,
            codon_table=_flat_codon_table(),
            skipped_codons_num=0,
        )
        assert result.detected_sites["motifs"] == 0
        assert set(result.detected_sites) == {"recombination", "slippage", "motifs"}


class TestEsoRevisionGuard:
    """widen_slippage_base_units declares num_base_units=2 for a repeat holding
    one widened unit, relying on the pinned ESO revision's modify_df_slippage
    emitting only the first of every two units. This guard is what makes an ESO
    upgrade that breaks the assumption say so."""

    def test_guard_is_silent_when_every_chunk_stays_inside(self):
        from modules.hotspot_avoidance.hotspot_avoidance_main import (
            _chunks_outside_detected_windows,
        )

        sequence = PLANTED_PREFIX + "CGCGCG" + PLANTED_SUFFIX
        widened, _ = widen_slippage_base_units(
            _slippage_df([30, 36, 2, "CGCGCG", 3, -4.2]), sequence
        )

        assert _chunks_outside_detected_windows(widened, {0: (30, 36)}) == []

    def test_guard_reports_a_chunk_that_escapes_the_detected_window(self):
        """Simulates the failure mode: if a future ESO emitted more chunks than
        the first of each pair, chunk 2 would land outside the repeat. Declaring
        4 units makes the CURRENT loop emit chunks 0 and 2, which reproduces
        exactly that shape without needing a modified ESO."""
        from modules.hotspot_avoidance.hotspot_avoidance_main import (
            _chunks_outside_detected_windows,
        )

        escaping = _slippage_df([30, 36, 6, "CGCGCG" * 4, 4, float("nan")])
        problems = _chunks_outside_detected_windows(escaping, {0: (30, 36)})

        assert len(problems) == 1
        assert "42-48" in problems[0]
        assert "pinned ESO revision" in problems[0]

    def test_guard_handles_an_empty_frame(self):
        from modules.hotspot_avoidance.hotspot_avoidance_main import (
            _chunks_outside_detected_windows,
        )

        assert _chunks_outside_detected_windows(_slippage_df(), {}) == []
