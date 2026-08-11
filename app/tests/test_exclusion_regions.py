import pytest

from modules.hotspot_avoidance.exclusion_regions import complement_regions
from modules.hotspot_avoidance.exclusion_regions import merge_regions
from modules.hotspot_avoidance.exclusion_regions import widen_to_codon_boundaries


class TestMergeRegions:
    def test_empty(self):
        assert merge_regions([]) == []

    def test_disjoint_regions_are_kept_apart_and_sorted(self):
        assert merge_regions([(50, 60), (10, 20)]) == [(10, 20), (50, 60)]

    def test_overlapping_regions_are_merged(self):
        assert merge_regions([(10, 30), (20, 40)]) == [(10, 40)]

    def test_touching_regions_are_merged(self):
        # (0, 45) and (45, 90) share no nucleotide, but as *editable* windows
        # they are contiguous, so collapsing them keeps the region list minimal.
        assert merge_regions([(0, 45), (45, 90)]) == [(0, 90)]

    def test_fully_contained_region_is_absorbed(self):
        assert merge_regions([(10, 100), (30, 40)]) == [(10, 100)]


class TestComplementRegions:
    def test_no_regions_means_everything_is_complement(self):
        assert complement_regions([], 100) == [(0, 100)]

    def test_single_middle_region(self):
        assert complement_regions([(30, 60)], 100) == [(0, 30), (60, 100)]

    def test_region_touching_the_start(self):
        assert complement_regions([(0, 30)], 100) == [(30, 100)]

    def test_region_touching_the_end(self):
        assert complement_regions([(60, 100)], 100) == [(0, 60)]

    def test_full_span_region_leaves_no_complement(self):
        assert complement_regions([(0, 100)], 100) == []

    def test_multiple_regions_are_merged_before_complementing(self):
        assert complement_regions([(20, 40), (30, 50)], 100) == [(0, 20), (50, 100)]

    def test_zero_length_sequence(self):
        assert complement_regions([], 0) == []


class TestWidenToCodonBoundaries:
    def test_already_aligned_region_is_unchanged(self):
        assert widen_to_codon_boundaries([(3, 9)]) == [(3, 9)]

    def test_start_rounds_down_and_end_rounds_up(self):
        assert widen_to_codon_boundaries([(4, 8)]) == [(3, 9)]

    def test_single_nucleotide_becomes_its_whole_codon(self):
        assert widen_to_codon_boundaries([(7, 8)]) == [(6, 9)]

    def test_multiple_regions_are_widened_independently(self):
        assert widen_to_codon_boundaries([(1, 2), (10, 11)]) == [(0, 3), (9, 12)]


def test_complement_of_widened_hotspots_locks_everything_else():
    # The property the whole locality guarantee rests on: a hotspot at (4, 8)
    # widens to (3, 9), and every other nucleotide of a 30nt sequence is locked.
    hotspots = widen_to_codon_boundaries([(4, 8)])
    assert complement_regions(hotspots, 30) == [(0, 3), (9, 30)]


import pandas as pd

from modules.hotspot_avoidance.exclusion_regions import build_exclusion_regions
from modules.hotspot_avoidance.exclusion_regions import hotspot_regions_from_detection


class TestHotspotRegionsFromDetection:
    def test_empty_detection_yields_no_regions(self):
        assert hotspot_regions_from_detection({}) == []

    def test_recombination_row_yields_both_members_of_the_pair(self):
        detection = {
            "df_recombination": pd.DataFrame(
                [{"sequence_1": "ACGTACGTACGTACGT", "start_1": 10, "end_1": 26,
                  "sequence_2": "ACGTACGTACGTACGT", "start_2": 100, "end_2": 116}]
            )
        }
        assert hotspot_regions_from_detection(detection) == [(10, 26), (100, 116)]

    def test_slippage_row_yields_one_region(self):
        detection = {
            "df_slippage": pd.DataFrame(
                [{"start": 30, "end": 42, "sequence": "AAAAAAAAAAAA", "length_base_unit": 1}]
            )
        }
        assert hotspot_regions_from_detection(detection) == [(30, 42)]

    def test_motif_end_index_is_inclusive_and_gets_plus_one(self):
        # eso.detection.methylation's end_index is the index of the motif's LAST
        # nucleotide, unlike every other detector. A GATC at 20..23 inclusive is
        # the exclusive-end region (20, 24). Getting this wrong is a real,
        # previously-shipped ESO bug (see eso/optimize.py's comment).
        detection = {
            "df_motifs": pd.DataFrame(
                [{"start_index": 20, "end_index": 23, "actual_site": "GATC"}]
            )
        }
        assert hotspot_regions_from_detection(detection) == [(20, 24)]

    def test_empty_dataframes_are_skipped(self):
        detection = {
            "df_recombination": pd.DataFrame(),
            "df_slippage": pd.DataFrame(),
            "df_motifs": pd.DataFrame(),
        }
        assert hotspot_regions_from_detection(detection) == []

    def test_all_three_types_are_combined_and_sorted(self):
        detection = {
            "df_recombination": pd.DataFrame(
                [{"sequence_1": "A" * 16, "start_1": 100, "end_1": 116,
                  "sequence_2": "A" * 16, "start_2": 200, "end_2": 216}]
            ),
            "df_slippage": pd.DataFrame(
                [{"start": 30, "end": 42, "sequence": "A" * 12, "length_base_unit": 1}]
            ),
            "df_motifs": pd.DataFrame(
                [{"start_index": 20, "end_index": 23, "actual_site": "GATC"}]
            ),
        }
        assert hotspot_regions_from_detection(detection) == [(20, 24), (30, 42), (100, 116), (200, 216)]


class TestBuildExclusionRegions:
    def test_no_hotspots_locks_the_entire_sequence(self):
        assert build_exclusion_regions([], sequence_length=90) == [(0, 90)]

    def test_single_hotspot_leaves_only_its_widened_window_editable(self):
        # (4, 8) widens to (3, 9); everything else is locked.
        assert build_exclusion_regions([(4, 8)], sequence_length=90) == [(0, 3), (9, 90)]

    def test_locked_prefix_is_never_editable(self):
        # build_exclusion_regions returns the LOCKED regions. The hotspot (0, 60)
        # survives only outside the locked prefix, so [45, 60) is the one editable
        # window - and everything on both sides of it stays locked.
        regions = build_exclusion_regions([(0, 60)], sequence_length=90, locked_prefix_length=45)
        assert regions == [(0, 45), (60, 90)]

    def test_hotspot_entirely_inside_the_locked_prefix_is_dropped(self):
        regions = build_exclusion_regions([(6, 12)], sequence_length=90, locked_prefix_length=45)
        assert regions == [(0, 90)]

    def test_adjacent_hotspots_produce_one_editable_window(self):
        # (3, 6) and (6, 9) are contiguous once widened, so the locked
        # complement must not contain a bogus zero-length gap between them.
        assert build_exclusion_regions([(3, 6), (6, 9)], sequence_length=30) == [(0, 3), (9, 30)]

    def test_overlapping_hotspots_are_merged(self):
        assert build_exclusion_regions([(3, 12), (9, 18)], sequence_length=30) == [(0, 3), (18, 30)]

    def test_hotspot_covering_everything_leaves_nothing_locked(self):
        assert build_exclusion_regions([(0, 30)], sequence_length=30) == []
