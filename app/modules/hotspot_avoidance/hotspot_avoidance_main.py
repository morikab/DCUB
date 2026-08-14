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
from eso.optimize import optimization_engine

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.hotspot_avoidance.dcub_score_adapter import build_dcub_codon_table
from modules.hotspot_avoidance.exclusion_regions import build_exclusion_regions
from modules.hotspot_avoidance.exclusion_regions import hotspot_regions_from_detection
from modules.shared_functions_and_vars import nt_to_aa
from modules.timer import Timer

logger = LoggerFactory.get_logger()
config = Configuration.get_config()

#: ESO's CustomScore unconditionally warns that it re-scores the whole ORF on
#: every trial mutation. That is a known, accepted cost of Approach B - not
#: something to show a biologist beside genuine "this hotspot could not be
#: cleared" warnings.
_ESO_PERFORMANCE_WARNING_FRAGMENT = "re-evaluates score_fn"


@dataclass
class HotspotPatchResult:
    sequence_before: str
    sequence_after: str
    num_edits: int
    detected_sites: typing.Dict[str, int]
    warnings: typing.List[str] = field(default_factory=list)

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "enabled": True,
            "sequence_before": self.sequence_before,
            "sequence_after": self.sequence_after,
            "num_edits": self.num_edits,
            "detected_sites": self.detected_sites,
            "warnings": self.warnings,
        }


def make_dcub_custom_score(
    codon_table: typing.Mapping[str, typing.Mapping[str, float]],
) -> typing.Callable[[str], float]:
    """Wrap a `{amino_acid: {codon: score}}` table as the `score(seq) -> float`
    callable ESO's CustomScore expects (higher is better).

    A codon with no table entry contributes 0, matching ESO's own
    "missing entry per codon -> 0" convention
    (single_codon_optimization_method._get_max_organism_attribute_value).
    """

    def dcub_custom_score(sequence: str) -> float:
        total = 0.0
        for index in range(0, len(sequence) - (len(sequence) % 3), 3):
            codon = sequence[index:index + 3].upper()
            amino_acid = nt_to_aa.get(codon)
            if amino_acid is None:
                continue
            total += codon_table.get(amino_acid, {}).get(codon, 0.0)
        return total

    return dcub_custom_score


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
        num_base_units = (end - start) // widened_base_unit

        if num_base_units < 2:
            # `modify_df_slippage` iterates range(0, num_base_units - 1, 2),
            # which is empty below 2, so such a row yields no constraints
            # anyway. Say so rather than dropping it silently.
            widening_warnings.append(
                f"Slippage site at {start}-{end} (base unit {length_base_unit}nt) is too short "
                f"to disrupt at codon resolution and was left unmodified."
            )
            continue

        widened_end = start + num_base_units * widened_base_unit
        widened_row = row.to_dict()
        widened_row.update({
            "length_base_unit": widened_base_unit,
            "num_base_units": num_base_units,
            "end": widened_end,
            # Recomputed, never carried over: the row must stay honest about
            # what is actually at these (possibly shortened) coordinates.
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


def patch_sequence(
    sequence: str,
    codon_table: typing.Mapping[str, typing.Mapping[str, float]],
    skipped_codons_num: int,
    compute_motifs: typing.Optional[bool] = None,
    common_motifs: typing.Optional[typing.List[str]] = None,
    recombination_mode: typing.Optional[str] = None,
    slippage_mode: typing.Optional[str] = None,
) -> HotspotPatchResult:
    """Detect hotspots in `sequence` and repair them without touching anything
    else. Returns the patched sequence plus everything the run summary reports.

    `compute_motifs`, `common_motifs`, `recombination_mode` and
    `slippage_mode` default to `app/modules/configuration.yaml`'s
    `HOTSPOT_AVOIDANCE` section when not supplied. The config is read here,
    inside the function body, rather than as a signature default - a
    signature default is evaluated once at import time and can never be
    monkeypatched afterwards.
    """
    hotspot_config = config["HOTSPOT_AVOIDANCE"]
    if compute_motifs is None:
        compute_motifs = hotspot_config["COMPUTE_MOTIFS"]
    if common_motifs is None:
        common_motifs = hotspot_config["COMMON_MOTIFS"]
    if recombination_mode is None:
        recombination_mode = hotspot_config["RECOMBINATION_MODE"]
    if slippage_mode is None:
        slippage_mode = hotspot_config["SLIPPAGE_MODE"]

    with Timer() as timer:
        extractor_kwargs = {
            "recombination_mode": recombination_mode,
            "slippage_mode": slippage_mode,
        }
        if compute_motifs:
            extractor_kwargs["common_motifs"] = list(common_motifs)
        detection = suspect_site_extractor(
            sequence,
            compute_motifs=compute_motifs,
            num_sites=np.inf,
            **extractor_kwargs,
        )

        detected_sites = {
            "recombination": int(len(detection.get("df_recombination", []))),
            "slippage": int(len(detection.get("df_slippage", []))),
            "motifs": int(len(detection.get("df_motifs", []))),
        }
        logger.info(f"Detected hypermutable sites: {detected_sites}")

        hotspot_regions = hotspot_regions_from_detection(detection)
        if not hotspot_regions:
            return HotspotPatchResult(
                sequence_before=sequence,
                sequence_after=sequence,
                num_edits=0,
                detected_sites=detected_sites,
                warnings=[],
            )

        exclusion_regions = build_exclusion_regions(
            hotspot_regions=hotspot_regions,
            sequence_length=len(sequence),
            locked_prefix_length=skipped_codons_num * 3,
        )

        # Detection counts above stay as DETECTED; widening only changes how the
        # same nucleotides are handed to DNAChisel for avoidance.
        df_slippage, slippage_warnings = widen_slippage_base_units(
            detection.get("df_slippage"),
            sequence,
        )

        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            final_sequence, _objectives_summary, num_edits = optimization_engine(
                sequence,
                # DCUB does not manage GC content, and with everything outside
                # the hotspots locked, an out-of-range window is unfixable -
                # ESO's retry loop would drop the constraint and emit a warning
                # that means nothing to the user. Disable it outright.
                mini_gc=0.0,
                maxi_gc=1.0,
                custom_score_fn=make_dcub_custom_score(codon_table),
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

        surfaced_warnings = slippage_warnings + [
            str(warning.message)
            for warning in caught_warnings
            if _ESO_PERFORMANCE_WARNING_FRAGMENT not in str(warning.message)
        ]
        for warning_message in surfaced_warnings:
            logger.warning(f"ESO hotspot avoidance: {warning_message}")

    logger.info(f"Hotspot avoidance made {num_edits} edits in {timer.elapsed_time}")
    return HotspotPatchResult(
        sequence_before=sequence,
        sequence_after=final_sequence,
        num_edits=int(num_edits),
        detected_sites=detected_sites,
        warnings=surfaced_warnings,
    )


class HotspotAvoidanceModule(object):
    @staticmethod
    def run_module(
        sequence: str,
        module_input: models.ModuleInput,
        optimization_cub_index: models.ORFOptimizationCubIndex,
        skipped_codons_num: int,
        compute_motifs: typing.Optional[bool] = None,
        common_motifs: typing.Optional[typing.List[str]] = None,
        recombination_mode: typing.Optional[str] = None,
        slippage_mode: typing.Optional[str] = None,
    ) -> HotspotPatchResult:
        codon_table = build_dcub_codon_table(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
        )
        return patch_sequence(
            sequence=sequence,
            codon_table=codon_table,
            skipped_codons_num=skipped_codons_num,
            compute_motifs=compute_motifs,
            common_motifs=common_motifs,
            recombination_mode=recombination_mode,
            slippage_mode=slippage_mode,
        )
