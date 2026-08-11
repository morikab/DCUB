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
from modules.ORF.zscore_optimization_method import _calculate_zscore_for_sequence
from modules.ORF.zscore_optimization_method import get_total_score
from modules.shared_functions_and_vars import change_all_codons_of_aa
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons

logger = LoggerFactory.get_logger()


def build_dcub_codon_table(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
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
        )

    if optimization_method.is_single_organism_optimization:
        return _table_from_wanted_organism_profile(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
        )

    if optimization_method.is_zscore_optimization:
        return _table_from_zscore(
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            sequence=sequence,
            skipped_codons_num=skipped_codons_num,
        )

    raise ValueError(
        f"Cannot build a codon preference table for optimization method {optimization_method}"
    )


def _table_from_codon_loss(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
) -> typing.Dict[str, typing.Dict[str, float]]:
    """single_codon_* methods pick the argMIN of a loss table. Negate it so the
    same codon becomes the argMAX, matching this module's higher-is-better
    contract.

    run_summary is deliberately None: ORFModule already wrote "orf_debug", and
    RunSummary raises KeyError on a duplicate key.
    """
    loss_table = _calculate_codons_loss(
        organisms=module_input.organisms,
        tuning_param=module_input.tuning_parameter,
        optimization_method=module_input.orf_optimization_method,
        optimization_cub_index=optimization_cub_index,
        run_summary=None,
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


def _table_from_zscore(
    module_input: models.ModuleInput,
    optimization_cub_index: models.ORFOptimizationCubIndex,
    sequence: str,
    skipped_codons_num: int,
) -> typing.Dict[str, typing.Dict[str, float]]:
    """Score every synonymous codon by the same computation
    optimize_sequence_by_zscore_bulk_aa runs each iteration - substitute the
    codon everywhere, z-score the resulting sequence, reduce to a total score -
    but run once, after the fact, against the FINAL candidate sequence.

    get_total_score is already higher-is-better, so no negation is needed.
    """
    optimization_method = module_input.orf_optimization_method

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

    if optimization_method.is_zscore_ratio_score_optimization:
        # _calculate_zscore_ratio_score takes a geometric mean, which is
        # undefined for the negative z-scores that standardization routinely
        # produces. Shift every candidate into a strictly positive range first,
        # exactly as optimize_sequence_by_zscore_bulk_aa does - without this the
        # whole table comes back nan.
        all_zscores = np.array(
            [point for zscore in codons_to_zscore.values() for point in zscore.all_scores]
        ).reshape(-1, 1)
        max_zscore_distance = pdist(all_zscores).max()
        min_zscore = all_zscores.min() - max_zscore_distance
        max_zscore = all_zscores.max() + max_zscore_distance
        for zscore in codons_to_zscore.values():
            zscore.normalize(min_zscore=min_zscore, max_zscore=max_zscore)

    table: typing.Dict[str, typing.Dict[str, float]] = {}
    for codon, zscore in codons_to_zscore.items():
        score = get_total_score(
            zscore=zscore,
            optimization_method=optimization_method,
            tuning_parameter=module_input.tuning_parameter,
        )
        table.setdefault(nt_to_aa[codon], {})[codon] = float(score)
    return table
