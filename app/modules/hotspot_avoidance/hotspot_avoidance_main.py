"""Run ESO's hypermutable-site detection on a DCUB-optimized sequence and
repair the detected hotspots LOCALLY, scoring replacements with DCUB's own
per-codon preference model.

The locality guarantee is mechanical rather than incidental: every nucleotide
outside a detected hotspot window is handed to DNAChisel as an `AvoidChanges`
constraint, so DCUB's differentially-optimized codon choices cannot drift no
matter which optimization method produced them.
"""

import math
import typing
import warnings
from dataclasses import dataclass
from dataclasses import field

import numpy as np
import pandas as pd
from eso import suspect_site_extractor
from eso.detection.slippage import modify_df_slippage
from eso.optimize import optimization_engine

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.hotspot_avoidance.dcub_score_adapter import build_dcub_score_fn
# Re-exported: this used to live here, and callers/tests import it from here.
from modules.hotspot_avoidance.dcub_score_adapter import make_dcub_custom_score  # noqa: F401
from modules.run_summary import RunSummary
from modules.hotspot_avoidance.exclusion_regions import build_exclusion_regions
from modules.hotspot_avoidance.exclusion_regions import hotspot_regions_from_detection
from modules.hotspot_avoidance.exclusion_regions import labeled_hotspot_regions_from_detection
from modules.shared_functions_and_vars import nt_to_aa
from modules.timer import Timer

logger = LoggerFactory.get_logger()
config = Configuration.get_config()

#: ESO's CustomScore unconditionally warns that it re-scores the whole ORF on
#: every trial mutation. That is a known, accepted cost of Approach B - not
#: something to show a biologist beside genuine "this hotspot could not be
#: cleared" warnings.
_ESO_PERFORMANCE_WARNING_FRAGMENT = "re-evaluates score_fn"

#: How many residual windows to name in the give-up warning before
#: summarising the rest as a count. A gene can end a run with dozens.
_MAX_RESIDUALS_LISTED = 5


@dataclass
class HotspotPatchResult:
    sequence_before: str
    sequence_after: str
    num_edits: int
    detected_sites: typing.Dict[str, int]
    warnings: typing.List[str] = field(default_factory=list)
    #: Every detected window as {"kind", "start", "end"}, 0-indexed with an
    #: exclusive end, in SEQUENCE_BEFORE coordinates. Repair is synonymous and
    #: therefore length-preserving, so these index sequence_after identically -
    #: but they mark what was DETECTED, which is not the same as what was
    #: edited: a window too narrow to disrupt is reported here and left alone.
    detected_regions: typing.List[typing.Dict[str, typing.Any]] = field(default_factory=list)
    #: How many repair rounds actually ran (0 when nothing was detected).
    rounds: int = 0
    #: Windows the verification pass still reported after the LAST round, in
    #: sequence_after coordinates. Empty means the final sequence came back
    #: clean; anything here is a site repair could not break, and is also
    #: surfaced in `warnings` so it reaches the results screen.
    residual_regions: typing.List[typing.Dict[str, typing.Any]] = field(default_factory=list)

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "enabled": True,
            "sequence_before": self.sequence_before,
            "sequence_after": self.sequence_after,
            "num_edits": self.num_edits,
            "detected_sites": self.detected_sites,
            "detected_regions": self.detected_regions,
            "rounds": self.rounds,
            "residual_regions": self.residual_regions,
            "warnings": self.warnings,
        }


def widen_slippage_base_units(
    df_slippage: typing.Optional[pd.DataFrame],
    sequence: str,
) -> typing.Tuple[typing.Optional[pd.DataFrame], typing.List[str]]:
    """Re-express SUB-CODON slippage base units - and only those, i.e.
    `length_base_unit` of 1 or 2 - as the smallest codon-width multiple, so
    they survive ESO's exclusion filtering.

    ESO's `exclusion_site_correcter` drops every avoidance row narrower than
    3nt (`df[df.start < df.end - 2]`) whenever exclusion regions are passed -
    and this module always passes them, since locking everything outside a
    hotspot IS the locality guarantee. `modify_df_slippage` emits rows exactly
    `length_base_unit` wide, so a homopolymer's 1nt rows are discarded before
    any constraint object exists, meaning not even a dropped-constraint
    warning is raised: the site is reported as detected and silently left
    unrepaired.

    Re-expressing a 15x"A" run as 5x"AAA" puts the same nucleotides under
    avoidance in rows wide enough to survive. Verified: 0 edits before, 2
    after, with flanks byte-identical.

    Base units of 3nt and wider are left exactly as detected. They already
    clear ESO's filter, so they never had the problem - and widening them
    would actively cause harm: it raises the repeat count needed for any
    repair from the detector's own minimum of 3 (see
    eso/detection/slippage.py `_generate_slippage_sites_current_subunit`) to
    6, and halves avoidance density above that. A 4nt-unit repeat like
    3x"ACGT" is repaired fine untouched, but becomes unrepairable if widened
    to 12nt units - `(42 - 30) // 12 == 1`, below the 2-unit minimum.

    Returns (widened_dataframe, warnings) - `warnings` names any run too
    short to disrupt at codon resolution, so it is reported rather than
    silently skipped.
    """
    if df_slippage is None or df_slippage.empty:
        return df_slippage, []

    rows = []
    widening_warnings = []
    for _, row in df_slippage.iterrows():
        length_base_unit = int(row["length_base_unit"])
        if length_base_unit >= 3:
            # Already wide enough to survive `df[df.start < df.end - 2]`.
            # See the docstring: widening these breaks working repairs.
            rows.append(row.to_dict())
            continue

        # 1 -> 3 and 2 -> 6: the smallest width that is both a whole number of
        # the detected base units (so the row still describes the real repeat)
        # and a whole number of codons.
        widened_base_unit = math.lcm(length_base_unit, 3)
        start = int(row["start"])
        end = int(row["end"])
        whole_units = (end - start) // widened_base_unit

        if whole_units < 1:
            # Not even one widened unit fits, so there is no chunk that both
            # covers whole codons AND lies inside the detected repeat. Nothing
            # correct to emit; say so rather than dropping it silently.
            widening_warnings.append(
                f"Slippage site at {start}-{end} (base unit {length_base_unit}nt) is too short "
                f"to disrupt at codon resolution and was left unmodified."
            )
            continue

        # ESO's modify_df_slippage iterates `range(0, num_base_units - 1, 2)`,
        # so it emits nothing at all below 2 units - which discarded every
        # repeat shorter than TWO widened units (12nt for a dinucleotide),
        # even though a single widened unit is already a perfectly good
        # avoidance target: it is codon-aligned, wider than ESO's 2nt filter,
        # and lies entirely inside the detected repeat.
        #
        # Declaring 2 units makes that loop run exactly once, emitting chunk 0
        # - the first widened_base_unit of the real repeat. Chunk 1 is never
        # emitted (the loop steps by 2), so the declared count being one higher
        # than the run really contains costs nothing.
        #
        # NOT the same as padding the window outward. A padded row would make
        # ESO avoid a pattern that straddles the repeat's flank, and DNAChisel
        # could then satisfy it by editing the FLANK and leaving the repeat
        # intact - measured on mCherry: padding one codon each side around
        # CGCGCG at 660-666 produced GAG->GAA, one reported edit, repeat fully
        # intact. Anchoring on the repeat produced GCG->GCC and destroyed it.
        # _assert_chunks_stay_inside below is what keeps that distinction
        # enforced rather than merely intended.
        num_base_units = max(2, whole_units)
        widened_end = start + whole_units * widened_base_unit
        widened_row = row.to_dict()
        widened_row.update({
            "length_base_unit": widened_base_unit,
            "num_base_units": num_base_units,
            "end": widened_end,
            # Recomputed, never carried over: the row must stay honest about
            # what is actually at these (possibly shortened) coordinates.
            # Sliced from the REAL repeat, never padded outward. Long enough
            # for modify_df_slippage to index chunk 0 out of it; when
            # num_base_units was raised to 2 this is still just the one real
            # widened unit, and chunk 1 is never requested.
            "sequence": sequence[start:widened_end],
            # Blanked for the same reason. ESO derives this from the base unit
            # and repeat count (-12.9 + 0.729n for 1nt units, -4.749 + 0.063n
            # otherwise), so carrying the detected value onto a re-parameterized
            # row leaves a triple that no longer satisfies ESO's own formula -
            # a trap for anyone who later reads it. Recomputing it from the
            # widened parameters would be worse: it would understate a real
            # homopolymer's risk by orders of magnitude (a 15x"A" run scores
            # -1.97, but 5x"AAA" through the non-homopolymer branch scores
            # -4.43) while looking authoritative. The widened row exists only
            # to generate AvoidPattern constraints; the true risk of the
            # physical site stays in the original detection dataframe, which
            # is what detected_sites is counted from.
            "log10_prob_slippage_ecoli": np.nan,
        })
        rows.append(widened_row)

    if not rows:
        # Preserve the column set - an empty DataFrame() would lose it, and
        # optimization_engine indexes these columns by name.
        return df_slippage.iloc[0:0].copy(), widening_warnings

    return pd.DataFrame(rows, columns=list(df_slippage.columns)), widening_warnings


def _chunks_outside_detected_windows(
    widened: pd.DataFrame,
    detected_windows: typing.Mapping[int, typing.Tuple[int, int]],
) -> typing.List[str]:
    """Verify every avoidance chunk ESO will derive from `widened` lies inside
    the window it came from.

    This is a version guard, not a sanity check. `widen_slippage_base_units`
    declares `num_base_units = 2` for a repeat that physically contains one
    widened unit, relying on ESO's
    `modify_df_slippage` iterating `range(0, num_base_units - 1, 2)` and so
    emitting only chunk 0. That holds for the ESO revision pinned in
    pyproject.toml. If a future revision changed the bound to
    `range(0, num_base_units, 2)`, chunk 1 would become an AvoidPattern over
    sequence that was never detected as a hotspot.

    The exclusion lock would contain the damage - that sequence is locked, so
    the constraint would be dropped as unsatisfiable - but it would be dropped
    noisily and for the wrong reason. Checking directly means an ESO upgrade
    that invalidates the assumption reports exactly that, instead of surfacing
    as confusing dropped-constraint warnings.

    Returns one message per out-of-bounds chunk; empty means the assumption
    still holds.
    """
    if widened.empty:
        return []

    problems = []
    for _, chunk in modify_df_slippage(widened.copy()).iterrows():
        chunk_start, chunk_end = int(chunk["start"]), int(chunk["end"])
        window = next(
            (
                (window_start, window_end)
                for window_start, window_end in detected_windows.values()
                if window_start <= chunk_start < window_end
            ),
            None,
        )
        if window is None or chunk_end > window[1]:
            problems.append(
                f"ESO derived an avoidance chunk at {chunk_start}-{chunk_end} that is not "
                f"contained in any detected slippage window. The pinned ESO revision's "
                f"modify_df_slippage emits only the first of every two base units; this "
                f"module relies on that. Re-check widen_slippage_base_units against the "
                f"ESO revision pinned in pyproject.toml."
            )
    return problems


def _repair_round(
    sequence: str,
    detection: typing.Mapping[str, typing.Any],
    score_fn: typing.Callable[[str], float],
    skipped_codons_num: int,
    padding_codons: int,
) -> typing.Tuple[str, typing.List[str]]:
    """One repair pass: lock everything outside `detection`'s windows, hand the
    windows to ESO as AvoidPattern constraints, return the patched sequence and
    whatever warnings it produced.

    `detection` must have been produced from `sequence` itself - the patterns
    ESO avoids are literal substrings at literal coordinates, so a detection
    taken from an earlier round's sequence would target strings that are no
    longer there.
    """
    exclusion_regions = build_exclusion_regions(
        hotspot_regions=hotspot_regions_from_detection(detection),
        sequence_length=len(sequence),
        locked_prefix_length=skipped_codons_num * 3,
        padding_codons=padding_codons,
    )

    # Detection counts stay as DETECTED; widening only changes how the same
    # nucleotides are handed to DNAChisel for avoidance.
    df_slippage, slippage_warnings = widen_slippage_base_units(
        detection.get("df_slippage"),
        sequence,
    )
    # Guard the ESO-revision assumption widen_slippage_base_units makes.
    # Cheap - the frame holds one row per detected site - and it turns an
    # ESO upgrade that breaks the assumption into a message naming the
    # cause, rather than a spray of dropped-constraint warnings.
    if df_slippage is not None and not df_slippage.empty:
        detected_windows = {
            index: (int(row["start"]), int(row["end"]))
            for index, row in enumerate(
                detection["df_slippage"].to_dict("records")
            )
        }
        slippage_warnings.extend(
            _chunks_outside_detected_windows(df_slippage, detected_windows)
        )

    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        final_sequence, _objectives_summary, _num_edits = optimization_engine(
            sequence,
            # DCUB does not manage GC content, and with everything outside
            # the hotspots locked, an out-of-range window is unfixable -
            # ESO's retry loop would drop the constraint and emit a warning
            # that means nothing to the user. Disable it outright.
            mini_gc=0.0,
            maxi_gc=1.0,
            custom_score_fn=score_fn,
            custom_score_minimize=False,
            df_recombination=detection.get("df_recombination"),
            df_slippage=df_slippage,
            df_motifs=detection.get("df_motifs"),
            # Pass the ORF explicitly: optimization_engine's own default is
            # ((len(seq) - 1) // 3) * 3, which drops the final codon of a
            # length-multiple-of-3 sequence (300 -> 297). DCUB sequences are
            # always whole codons.
            orf_regions=[(0, len(sequence))],
            exclusion_regions=exclusion_regions,
        )

    round_warnings = slippage_warnings + [
        str(warning.message)
        for warning in caught_warnings
        if _ESO_PERFORMANCE_WARNING_FRAGMENT not in str(warning.message)
    ]
    return final_sequence, round_warnings


def _describe_residual_regions(
    residual_regions: typing.List[typing.Dict[str, typing.Any]],
    rounds: int,
) -> str:
    """The warning shown when repair ran out of rounds with sites still
    standing. Positions are 1-indexed and inclusive, matching how the results
    screen lists them and how a biologist reads a position."""
    listed = ", ".join(
        f"{region['kind']} {region['start'] + 1}-{region['end']}"
        for region in residual_regions[:_MAX_RESIDUALS_LISTED]
    )
    if len(residual_regions) > _MAX_RESIDUALS_LISTED:
        listed += f", and {len(residual_regions) - _MAX_RESIDUALS_LISTED} more"
    return (
        f"{len(residual_regions)} hypermutable site(s) still detected after "
        f"{rounds} repair round(s): {listed}. Raise "
        f"HOTSPOT_AVOIDANCE.MAX_REPAIR_ROUNDS in app/modules/configuration.yaml "
        f"to allow more attempts."
    )


def patch_sequence(
    sequence: str,
    skipped_codons_num: int,
    codon_table: typing.Optional[typing.Mapping[str, typing.Mapping[str, float]]] = None,
    score_fn: typing.Optional[typing.Callable[[str], float]] = None,
    compute_motifs: typing.Optional[bool] = None,
    common_motifs: typing.Optional[typing.List[str]] = None,
    recombination_mode: typing.Optional[str] = None,
    slippage_mode: typing.Optional[str] = None,
    max_rounds: typing.Optional[int] = None,
) -> HotspotPatchResult:
    """Detect hotspots in `sequence` and repair them without touching anything
    else. Returns the patched sequence plus everything the run summary reports.

    Repair is ITERATIVE: every round is followed by a fresh detection pass on
    the round's own output, and another round runs if anything is still
    standing, up to `max_rounds` attempts. Neither ESO nor a single round
    verifies its own work - ESO's AvoidPattern constraints guarantee the exact
    detected string is gone from the exact detected window, which is not the
    same as the window having stopped being hypermutable. Each retry also
    unlocks one more codon on each side of every window (see
    `build_exclusion_regions`), because the commonest survivor is a repeat
    whose remaining copy sits just outside the previous round's reach.

    `compute_motifs`, `common_motifs`, `recombination_mode`, `slippage_mode`
    and `max_rounds` default to `app/modules/configuration.yaml`'s
    `HOTSPOT_AVOIDANCE` section when not supplied. The config is read here,
    inside the function body, rather than as a signature default - a
    signature default is evaluated once at import time and can never be
    monkeypatched afterwards.

    Supply exactly one of `codon_table` or `score_fn`. A table is the natural
    form for the two per-codon method families; `score_fn` exists because the
    zscore family has no per-codon decomposition - its score is a property of
    the whole sequence - and is scored exactly instead of approximated.
    """
    if (codon_table is None) == (score_fn is None):
        raise ValueError("patch_sequence needs exactly one of codon_table or score_fn")
    if score_fn is None:
        score_fn = make_dcub_custom_score(codon_table)

    hotspot_config = config["HOTSPOT_AVOIDANCE"]
    if compute_motifs is None:
        compute_motifs = hotspot_config["COMPUTE_MOTIFS"]
    if common_motifs is None:
        common_motifs = hotspot_config["COMMON_MOTIFS"]
    if max_rounds is None:
        max_rounds = hotspot_config["MAX_REPAIR_ROUNDS"]
    # A round count below 1 would skip repair altogether while still reporting
    # detected sites - silently doing nothing is the one outcome worth ruling
    # out here, so floor it at a single pass.
    max_rounds = max(1, int(max_rounds))

    with Timer() as timer:
        # Detector modes are omitted unless explicitly overridden, so ESO's own
        # defaults apply rather than values restated here. Naming them would
        # pin whatever ESO's defaults happened to be when this was written.
        extractor_kwargs = {}
        if recombination_mode is not None:
            extractor_kwargs["recombination_mode"] = recombination_mode
        if slippage_mode is not None:
            extractor_kwargs["slippage_mode"] = slippage_mode

        if compute_motifs:
            # COMMON_MOTIFS is NOT redundant with ESO's default, which is None.
            # suspect_site_extractor guards with `if common_motifs:`, so a None
            # default with no motifs_path leaves the motif list empty and
            # find_motif_sites returns an empty frame - motif detection would be
            # silently enabled and do nothing. Verified directly: with
            # compute_motifs=True and common_motifs omitted, a sequence carrying
            # 5x GATC and 4x CCAGG reports 0 motif rows; with ["dam", "dcm"] it
            # reports 12.
            if not common_motifs:
                raise ValueError(
                    "Motif detection is enabled but no motifs are configured. Set "
                    "HOTSPOT_AVOIDANCE.COMMON_MOTIFS in app/modules/configuration.yaml "
                    "(e.g. [\"dam\", \"dcm\"]), or set COMPUTE_MOTIFS to False."
                )
            extractor_kwargs["common_motifs"] = list(common_motifs)

        def detect(target_sequence: str) -> typing.Dict[str, typing.Any]:
            return suspect_site_extractor(
                target_sequence,
                compute_motifs=compute_motifs,
                num_sites=np.inf,
                **extractor_kwargs,
            )

        detection = detect(sequence)
        detected_sites = {
            "recombination": int(len(detection.get("df_recombination", []))),
            "slippage": int(len(detection.get("df_slippage", []))),
            "motifs": int(len(detection.get("df_motifs", []))),
        }
        logger.info(f"Detected hypermutable sites: {detected_sites}")

        # What the FIRST detection found, in sequence_before coordinates - this
        # is what the results screen highlights. Sites that only appear in a
        # later round are consequences of repair, not properties of the
        # incoming sequence, and are reported as residuals instead.
        detected_regions = labeled_hotspot_regions_from_detection(detection)
        if not detected_regions:
            return HotspotPatchResult(
                sequence_before=sequence,
                sequence_after=sequence,
                num_edits=0,
                detected_sites=detected_sites,
                warnings=[],
                detected_regions=detected_regions,
                rounds=0,
                residual_regions=[],
            )

        current_sequence = sequence
        collected_warnings: typing.List[str] = []
        residual_regions = detected_regions
        rounds = 0
        for round_index in range(max_rounds):
            current_sequence, round_warnings = _repair_round(
                sequence=current_sequence,
                detection=detection,
                score_fn=score_fn,
                skipped_codons_num=skipped_codons_num,
                # Round 0 keeps the tight, hotspot-only window; each retry
                # unlocks one more codon on each side.
                padding_codons=round_index,
            )
            collected_warnings.extend(round_warnings)
            rounds += 1

            # The verification pass. It runs even on the last round, so a
            # sequence that ships with sites still in it says so.
            detection = detect(current_sequence)
            residual_regions = labeled_hotspot_regions_from_detection(detection)
            if not residual_regions:
                break
            logger.info(
                f"Repair round {rounds} left {len(residual_regions)} site(s) standing; "
                f"retrying with {round_index + 1} codon(s) of padding"
                if rounds < max_rounds
                else f"Repair round {rounds} left {len(residual_regions)} site(s) standing"
            )

        if residual_regions:
            collected_warnings.append(
                _describe_residual_regions(residual_regions, rounds)
            )

        # Rounds repeat the same constraints, so the same message can arrive
        # several times; keep first occurrences only, in order.
        surfaced_warnings = list(dict.fromkeys(collected_warnings))
        for warning_message in surfaced_warnings:
            logger.warning(f"ESO hotspot avoidance: {warning_message}")

    # Counted against the ORIGINAL sequence rather than summed over rounds:
    # repair is synonymous and length-preserving, so this is an exact
    # position-wise diff, and a position edited in two rounds counts once.
    num_edits = sum(
        1 for before, after in zip(sequence, current_sequence) if before != after
    )
    logger.info(
        f"Hotspot avoidance made {num_edits} edits over {rounds} round(s) "
        f"in {timer.elapsed_time}"
    )
    return HotspotPatchResult(
        sequence_before=sequence,
        sequence_after=current_sequence,
        num_edits=num_edits,
        detected_sites=detected_sites,
        warnings=surfaced_warnings,
        detected_regions=detected_regions,
        rounds=rounds,
        residual_regions=residual_regions,
    )


class HotspotAvoidanceModule(object):
    @staticmethod
    def run_module(
        sequence: str,
        module_input: models.ModuleInput,
        optimization_cub_index: models.ORFOptimizationCubIndex,
        skipped_codons_num: int,
        run_summary: RunSummary,
        compute_motifs: typing.Optional[bool] = None,
        common_motifs: typing.Optional[typing.List[str]] = None,
        recombination_mode: typing.Optional[str] = None,
        slippage_mode: typing.Optional[str] = None,
        max_rounds: typing.Optional[int] = None,
    ) -> HotspotPatchResult:
        score_fn = build_dcub_score_fn(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
            run_summary=run_summary,
        )
        return patch_sequence(
            sequence=sequence,
            score_fn=score_fn,
            skipped_codons_num=skipped_codons_num,
            compute_motifs=compute_motifs,
            common_motifs=common_motifs,
            recombination_mode=recombination_mode,
            slippage_mode=slippage_mode,
            max_rounds=max_rounds,
        )
