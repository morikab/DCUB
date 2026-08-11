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
