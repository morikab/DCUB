"""Expose DCUB's per-codon preference model as a plain
`{amino_acid: {codon: score}}` table, normalized so HIGHER IS ALWAYS BETTER.

This is the boundary that lets the hotspot-avoidance module hand ESO a scoring
function reflecting DCUB's wanted/unwanted-organism tradeoffs, instead of
ESO's generic CAI/tAI tables - so a forced substitution inside a hotspot
window still respects the differential optimization DCUB just performed.

Normalizing the direction here means nothing downstream needs to know which
DCUB method family produced the table.
"""

import typing

import numpy as np
from scipy.spatial.distance import pdist

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.ORF.single_codon_optimization_method import _calculate_codons_loss
from modules.ORF.single_codon_optimization_method import orf_detailed_summary_key
from modules.configuration import Configuration
from modules.ORF.zscore_optimization_method import _calculate_zscore_for_sequence
from modules.ORF.zscore_optimization_method import get_total_score
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import change_all_codons_of_aa
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons

logger = LoggerFactory.get_logger()
config = Configuration.get_config()

#: Subtracted once per rare codon appearing in a scored sequence. It has to
#: outweigh any achievable difference in the underlying score so that a codon
#: clearing the frequency floor always beats one that does not - which is
#: exactly how DCUB's own _get_optimal_codon treats the floor: a hard skip, not
#: a tiebreak. Safe to make this large: objectives are SOFT in DNAChisel
#: (resolve_constraints runs first and to completion), so no penalty can ever
#: prevent a hotspot from being repaired - it only decides WHICH legal codon
#: gets used.
RARE_CODON_PENALTY = 1.0e6


def rare_codons(organisms: typing.Sequence[models.Organism]) -> typing.FrozenSet[str]:
    """Codons DCUB's own `_get_optimal_codon` would refuse to select: those
    whose mean frequency across the WANTED organisms falls below
    `config["ORF"]["FREQUENCY_OPTIMIZATION_THRESHOLD"]`.

    Amino acids whose every synonymous codon is below the floor are excluded,
    mirroring `_get_optimal_codon`'s fallback ("could not find codon that
    satisfies minimal average frequency ... using the codon with the minimal
    loss score"). Penalizing those would be pointless - there is no compliant
    alternative to steer toward - and would distort scores between amino acids
    for no benefit.
    """
    wanted = [organism for organism in organisms if organism.is_optimized]
    if not wanted:
        return frozenset()

    threshold = config["ORF"]["FREQUENCY_OPTIMIZATION_THRESHOLD"]

    def below_floor(codon: str) -> bool:
        frequencies = [organism.codon_frequencies.get(codon, 0.0) for organism in wanted]
        return (sum(frequencies) / len(frequencies)) < threshold

    penalized = set()
    for codons in synonymous_codons.values():
        offending = {codon for codon in codons if below_floor(codon)}
        if len(offending) < len(codons):
            penalized |= offending
    return frozenset(penalized)


def _rare_codon_count(sequence: str, penalized: typing.FrozenSet[str]) -> int:
    if not penalized:
        return 0
    return sum(
        1
        for index in range(0, len(sequence) - (len(sequence) % 3), 3)
        if sequence[index:index + 3].upper() in penalized
    )


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


def build_dcub_codon_table(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
    run_summary: RunSummary,
) -> typing.Dict[str, typing.Dict[str, float]]:
    """Return DCUB's per-codon preferences as a higher-is-better table.

    `sequence` and `skipped_codons_num` are only used by the zscore family,
    whose scores are a property of the whole sequence rather than of a codon
    in isolation.
    """
    optimization_method = module_input.orf_optimization_method

    if optimization_method.is_single_codon_optimization:
        return _table_from_codon_loss(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            run_summary=run_summary,
        )

    if optimization_method.is_single_organism_optimization:
        return _table_from_wanted_organism_profile(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
        )

    raise ValueError(
        f"Cannot build a codon preference table for optimization method "
        f"{optimization_method}. The zscore family has no per-codon table - its "
        f"score is a property of the whole sequence; use build_dcub_score_fn."
    )


def _table_from_codon_loss(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    run_summary: RunSummary,
) -> typing.Dict[str, typing.Dict[str, float]]:
    """single_codon_* methods pick the argMIN of a loss table. Negate it so the
    same codon becomes the argMAX, matching this module's higher-is-better
    contract.

    Recorded under its own "hotspot_avoidance_detailed_<index>" key rather
    than overwriting the ORF module's. The two tables should be identical -
    same organisms, tuning parameter, method and CUB index - so keeping them
    as separate entries makes any divergence between the two stages visible
    in the run summary instead of silently replacing one with the other.
    """
    loss_table = _calculate_codons_loss(
        organisms=module_input.organisms,
        tuning_param=module_input.tuning_parameter,
        optimization_method=module_input.orf_optimization_method,
        optimization_cub_index=optimization_cub_index,
        run_summary=run_summary,
        summary_key=orf_detailed_summary_key(
            optimization_cub_index, stage="hotspot_avoidance"
        ),
    )
    return {
        amino_acid: {codon: -loss for codon, loss in codon_losses.items()}
        for amino_acid, codon_losses in loss_table.items()
    }


def _table_from_wanted_organism_profile(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
) -> typing.Dict[str, typing.Dict[str, float]]:
    """single_wanted_organism picks the argMAX of the wanted organism's CUB
    profile directly (see ORF/single_organism_optimization_method._get_optimal_codons),
    so that profile IS the higher-is-better table.

    The per-amino-acid dict is built by walking `profile.items()` (falling
    back to `synonymous_codons` order only for codons the profile omits)
    rather than walking `synonymous_codons` order directly, so that ties are
    broken identically to `_get_optimal_codons`'s own `candidate_codons`
    construction - `_get_optimal_codons` builds its candidate dict from
    `cub_profile.items()` too. Building from `synonymous_codons` order
    instead makes `max()` resolve ties (equal-weight codons) differently
    from `_get_optimal_codons`, so the two argmax picks silently diverge
    whenever a profile has tied weights - verified against this module's own
    test fixture, which ties every amino acid but Cysteine at 0.2 and hits
    this exact divergence for 14 of 21 amino acids without this fix.
    """
    wanted_organisms = [organism for organism in module_input.organisms if organism.is_optimized]
    if len(wanted_organisms) != 1:
        logger.warning(
            f"Number of wanted organisms is {len(wanted_organisms)} for optimization method "
            f"{module_input.orf_optimization_method}. Building the codon table from the first "
            f"wanted organism, matching _get_optimal_codons."
        )
    wanted_organism = wanted_organisms[0]
    profile = getattr(wanted_organism, f"{optimization_cub_index.value.lower()}_profile", {}) or {}

    table = {}
    for amino_acid, codons in synonymous_codons.items():
        codon_scores = {codon: float(weight) for codon, weight in profile.items() if codon in codons}
        for codon in codons:
            codon_scores.setdefault(codon, 0.0)
        table[amino_acid] = codon_scores
    return table


def _zscore_per_codon_sweep(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
) -> typing.Dict[str, typing.Any]:
    """z-score every "what if EVERY codon of this amino acid became X" variant
    of `sequence`, the same sweep optimize_sequence_by_zscore_bulk_aa performs
    each iteration.

    Only used to derive stable normalization bounds for the ratio family - the
    score itself is computed exactly, per trial sequence, by
    _exact_zscore_score_fn.
    """
    codons_to_zscore = {}
    for codon in nt_to_aa:
        candidate_sequence = change_all_codons_of_aa(
            seq=sequence,
            selected_codon=codon,
            skipped_codons_num=skipped_codons_num,
        )
        codons_to_zscore[codon] = _calculate_zscore_for_sequence(
            sequence=candidate_sequence,
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            skipped_codons_num=skipped_codons_num,
        )
    return codons_to_zscore


def _ratio_normalization_bounds(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
) -> typing.Tuple[float, float]:
    """Fixed (min, max) bounds for shifting z-scores into a strictly positive
    range before the ratio family's geometric mean.

    _calculate_zscore_ratio_score takes a geometric mean, undefined for the
    negative z-scores standardization routinely produces, so without this every
    score comes back nan. The bounds are derived ONCE, from the per-codon sweep
    over the candidate sequence, and then held fixed for every trial evaluation
    - the same reason optimize_sequence_by_zscore_bulk_aa carries running
    bounds across its convergence loop: a normalization that moves with each
    trial makes successive scores incomparable, which is exactly what an
    optimizer must be able to compare.
    """
    codons_to_zscore = _zscore_per_codon_sweep(
        module_input=module_input,
        optimization_cub_index=optimization_cub_index,
        sequence=sequence,
        skipped_codons_num=skipped_codons_num,
    )
    all_zscores = np.array(
        [point for zscore in codons_to_zscore.values() for point in zscore.all_scores]
    ).reshape(-1, 1)
    max_zscore_distance = pdist(all_zscores).max()
    return (
        float(all_zscores.min() - max_zscore_distance),
        float(all_zscores.max() + max_zscore_distance),
    )


def _exact_zscore_score_fn(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
) -> typing.Callable[[str], float]:
    """DCUB's REAL zscore objective, evaluated on each trial sequence.

    This replaces a per-codon approximation that asked "what if every codon of
    this amino acid became X?" and scored that against the final candidate,
    once. DCUB's actual zscore optimization is positional and iterative, so the
    proxy table computed a different, structurally unrelated quantity and then
    ranked codons by it.

    Exact evaluation is affordable precisely because of the locality lock: only
    codons inside a detected hotspot are editable, so DNAChisel's objective
    pass runs few trials. Measured at 0.383 ms per call on a 711nt gene - and
    the discarded proxy table already spent 64 of these calls building itself,
    so this is cheaper than what it replaces in any run with fewer than 64
    trials.
    """
    optimization_method = module_input.orf_optimization_method
    bounds = None
    if optimization_method.is_zscore_ratio_score_optimization:
        bounds = _ratio_normalization_bounds(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
        )

    def exact_zscore_score(candidate: str) -> float:
        zscore = _calculate_zscore_for_sequence(
            sequence=candidate,
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            skipped_codons_num=skipped_codons_num,
        )
        if bounds is not None:
            zscore.normalize(min_zscore=bounds[0], max_zscore=bounds[1])
        return float(get_total_score(
            zscore=zscore,
            optimization_method=optimization_method,
            tuning_parameter=module_input.tuning_parameter,
        ))

    return exact_zscore_score


def build_dcub_score_fn(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
    run_summary: RunSummary,
) -> typing.Callable[[str], float]:
    """DCUB's preference model as the `score(seq) -> float` callable ESO's
    CustomScore expects (higher is better), for ANY optimization method.

    Two families reduce to a per-codon table; the zscore family is scored
    exactly, on the trial sequence. Both then get the same rare-codon penalty
    applied on top, so the frequency floor DCUB's own `_get_optimal_codon`
    enforces holds no matter which method produced the base score. Without it
    the optimizer could install a codon DCUB itself refused to select for being
    too rare in the wanted hosts - verified: with Cys weighted so TGT has the
    minimal loss but only 0.02 mean frequency in the wanted host, DCUB picks
    TGC while an unpenalized table picks TGT.
    """
    optimization_method = module_input.orf_optimization_method

    if optimization_method.is_zscore_optimization:
        base_score_fn = _exact_zscore_score_fn(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
        )
    else:
        codon_table = build_dcub_codon_table(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
            run_summary=run_summary,
        )
        base_score_fn = make_dcub_custom_score(codon_table)

    return apply_rare_codon_penalty(base_score_fn, module_input.organisms)


def apply_rare_codon_penalty(
    score_fn: typing.Callable[[str], float],
    organisms: typing.Sequence[models.Organism],
) -> typing.Callable[[str], float]:
    """Wrap `score_fn` so each rare codon in the scored sequence costs
    RARE_CODON_PENALTY.

    Applied as a term on the score rather than by removing codons from a table
    because it has to work for the zscore family too, which has no table. The
    penalty is per-codon, and DNAChisel edits one codon at a time, so a
    compliant synonymous alternative always outranks a rare one by at least
    RARE_CODON_PENALTY.
    """
    penalized = rare_codons(organisms)
    if not penalized:
        return score_fn

    logger.info(
        f"{len(penalized)} codon(s) fall below the wanted-host frequency floor "
        f"and will be avoided where a synonymous alternative exists: "
        f"{sorted(penalized)}"
    )

    def penalized_score(candidate: str) -> float:
        return score_fn(candidate) - RARE_CODON_PENALTY * _rare_codon_count(candidate, penalized)

    return penalized_score
