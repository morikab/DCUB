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
