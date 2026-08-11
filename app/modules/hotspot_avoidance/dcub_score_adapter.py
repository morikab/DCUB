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

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.ORF.single_codon_optimization_method import _calculate_codons_loss
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
