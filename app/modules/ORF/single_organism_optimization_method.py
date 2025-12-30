import typing

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules import shared_functions_and_vars
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.timer import Timer

logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
def optimize_sequence(
        target_gene: str,
        organisms: typing.Sequence[models.Organism],
        optimization_method: models.ORFOptimizationMethod,
        optimization_cub_index: models.ORFOptimizationCubIndex,
        skipped_codons_num: int,
        run_summary: RunSummary,
) -> str:
    with Timer() as timer:
        aa_to_optimal_codon = _get_optimal_codons(
            organisms=organisms,
            optimization_method=optimization_method,
            optimization_cub_index=optimization_cub_index,
        )

        target_protein = shared_functions_and_vars.translate(target_gene)

        skipped_codons_size_in_nt = skipped_codons_num * 3
        optimized_sequence = target_gene[:skipped_codons_size_in_nt]
        for i in range(skipped_codons_num, len(target_protein)):
            aa = target_protein[i]
            optimal_codon = aa_to_optimal_codon[aa]
            optimized_sequence += optimal_codon

        if target_protein.endswith("_") and optimization_cub_index.is_trna_adaptation_index:
            # There is no point in optimizing stop codon by tAI, so leaving the original codon
            optimized_sequence = optimized_sequence[:-3]
            optimized_sequence += target_gene[-3:]

    orf_summary = {
        "orf_module_input_sequence": target_gene,
        "optimized_sequence": optimized_sequence,
        "aa_to_optimal_codon": aa_to_optimal_codon,
        "run_time": timer.elapsed_time,
    }
    run_summary.add_to_run_summary("orf", orf_summary)
    return optimized_sequence

# --------------------------------------------------------------
def _get_optimal_codons(
    organisms: typing.Sequence[models.Organism],
    optimization_method: models.ORFOptimizationMethod,
    optimization_cub_index: models.ORFOptimizationCubIndex,
) -> dict[str, str]:
    wanted_organisms = [o for o in organisms if o.is_optimized]
    wanted_organisms_count = len(wanted_organisms)
    if wanted_organisms_count != 1:
        logger.warning(f"Number of wanted organisms is {wanted_organisms_count} for the current optimization method "
                       f"{optimization_method}. Optimizing based on the first wanted organism.")
    wanted_organism = wanted_organisms[0]
    aa_to_optimal_codon = {}
    for aa, codons in shared_functions_and_vars.synonymous_codons.items():
        cub_profile = getattr(wanted_organism, f"{optimization_cub_index.value.lower()}_profile", {})
        candidate_codons = {codon: weight for codon, weight in cub_profile.items() if codon in codons}
        optimal_codon = max(candidate_codons, key=candidate_codons.get)
        aa_to_optimal_codon[aa] = optimal_codon

    return aa_to_optimal_codon
