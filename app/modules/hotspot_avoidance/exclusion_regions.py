"""Interval algebra for turning ESO's detected hotspot windows into the
`exclusion_regions` (locked, unmodifiable) list that `eso.optimize` consumes.

Every region here is 0-indexed with an EXCLUSIVE end - `seq[start:end]` -
matching ESO's own convention. Two touching regions such as (0, 45) and
(45, 90) share no actual nucleotide.
"""

import typing


def merge_regions(regions: typing.Iterable[typing.Tuple[int, int]]) -> typing.List[typing.Tuple[int, int]]:
    """Sort and coalesce overlapping or touching regions.

    Touching regions ARE merged here, unlike in overlap *detection*: two
    editable windows that run (0, 45) and (45, 90) are contiguously editable,
    so emitting them as one (0, 90) window keeps the constraint list minimal
    without changing which nucleotides are editable.
    """
    ordered = sorted((int(start), int(end)) for start, end in regions if int(end) > int(start))
    if not ordered:
        return []

    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def complement_regions(
    regions: typing.Iterable[typing.Tuple[int, int]],
    sequence_length: int,
) -> typing.List[typing.Tuple[int, int]]:
    """Return everything in [0, sequence_length) that `regions` does not cover."""
    complement = []
    cursor = 0
    for start, end in merge_regions(regions):
        start = max(0, min(start, sequence_length))
        end = max(0, min(end, sequence_length))
        if start > cursor:
            complement.append((cursor, start))
        cursor = max(cursor, end)

    if cursor < sequence_length:
        complement.append((cursor, sequence_length))
    return complement


def widen_to_codon_boundaries(
    regions: typing.Iterable[typing.Tuple[int, int]],
) -> typing.List[typing.Tuple[int, int]]:
    """Expand each region outward to whole codons.

    A hotspot rarely lands on a codon boundary, but the only edit ESO can
    legally make is a synonymous codon swap. Widening outward means the
    editable window always contains whole codons, so there is genuinely
    something to change; widening *inward* could produce a window with no
    complete codon in it at all.
    """
    return [((int(start) // 3) * 3, -(-int(end) // 3) * 3) for start, end in regions]


def labeled_hotspot_regions_from_detection(
    detection: typing.Mapping[str, typing.Any],
) -> typing.List[typing.Dict[str, typing.Any]]:
    """Flatten `eso.suspect_site_extractor`'s dict of dataframes into a sorted
    list of `{"kind", "start", "end"}` hotspot windows.

    Coordinates are 0-indexed with an EXCLUSIVE end, into the sequence
    detection ran on - i.e. the pre-repair sequence.

    Each detector reports its coordinates differently:
      - df_recombination: a PAIR of regions per row (start_1/end_1 and
        start_2/end_2) - both are hotspot windows, since breaking the
        near-duplicate relationship may require editing either one.
      - df_slippage: one region per row, start/end, exclusive end.
      - df_motifs: start_index/end_index, where end_index is INCLUSIVE - the
        index of the motif's last nucleotide. The +1 below converts it to this
        codebase's exclusive-end convention. Omitting it makes every motif
        window one nucleotide too short to contain its own motif, which in ESO
        itself once made motif avoidance silently do nothing.

    `kind` is one of "recombination" / "slippage" / "motifs", matching the keys
    of `HotspotPatchResult.detected_sites`.
    """
    regions: typing.List[typing.Dict[str, typing.Any]] = []

    df_recombination = detection.get("df_recombination")
    if df_recombination is not None and not df_recombination.empty:
        for _, row in df_recombination.iterrows():
            regions.append({"kind": "recombination",
                            "start": int(row["start_1"]), "end": int(row["end_1"])})
            regions.append({"kind": "recombination",
                            "start": int(row["start_2"]), "end": int(row["end_2"])})

    df_slippage = detection.get("df_slippage")
    if df_slippage is not None and not df_slippage.empty:
        for _, row in df_slippage.iterrows():
            regions.append({"kind": "slippage",
                            "start": int(row["start"]), "end": int(row["end"])})

    df_motifs = detection.get("df_motifs")
    if df_motifs is not None and not df_motifs.empty:
        for _, row in df_motifs.iterrows():
            regions.append({"kind": "motifs",
                            "start": int(row["start_index"]),
                            "end": int(row["end_index"]) + 1})

    return sorted(regions, key=lambda region: (region["start"], region["end"], region["kind"]))


def hotspot_regions_from_detection(
    detection: typing.Mapping[str, typing.Any],
) -> typing.List[typing.Tuple[int, int]]:
    """The same windows as `labeled_hotspot_regions_from_detection`, reduced to
    sorted (start, end) tuples - what `build_exclusion_regions` consumes.
    """
    return sorted(
        (region["start"], region["end"])
        for region in labeled_hotspot_regions_from_detection(detection)
    )


def build_exclusion_regions(
    hotspot_regions: typing.Iterable[typing.Tuple[int, int]],
    sequence_length: int,
    locked_prefix_length: int = 0,
) -> typing.List[typing.Tuple[int, int]]:
    """Return the regions ESO must NOT modify: everything except the
    codon-boundary-widened hotspot windows, plus the initiation-optimized
    prefix.

    This is the mechanical locality guarantee. Rather than trusting the
    objective not to wander, every nucleotide outside a detected hotspot is
    handed to DNAChisel as an `AvoidChanges` constraint, so DCUB's chosen
    codons cannot drift no matter which optimization method produced them.
    """
    editable = merge_regions(widen_to_codon_boundaries(hotspot_regions))

    # Clip the editable windows out of the initiation-optimized prefix. Doing
    # it here (rather than adding the prefix back as an exclusion afterwards)
    # means the complement below produces one contiguous leading locked region
    # instead of two abutting ones.
    if locked_prefix_length > 0:
        editable = [
            (max(start, locked_prefix_length), end)
            for start, end in editable
            if end > locked_prefix_length
        ]

    return complement_regions(editable, sequence_length)
