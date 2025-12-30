import typing
from functools import partial

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import synonymous_codon_permutation
from modules.shared_functions_and_vars import validate_module_output

from .single_codon_optimization_method import optimize_sequence as optimize_sequence_by_single_codon
from .single_organism_optimization_method import optimize_sequence as optimize_sequence_for_single_organism
from .zscore_optimization_method import optimize_sequence_by_zscore_single_aa
from .zscore_optimization_method import optimize_sequence_by_zscore_bulk_aa

config = Configuration.get_config()
logger = LoggerFactory.get_logger()


class ORFModule(object):
    @staticmethod
    def run_module(
            module_input: models.ModuleInput,
            optimization_cub_index: models.ORFOptimizationCubIndex,
            optimization_method: models.ORFOptimizationMethod,
            skipped_codons_num: int,
            run_summary: RunSummary,
    ) -> typing.Sequence[str]:
        logger.info("\n##########################")
        logger.info("# ORF #")
        logger.info("##########################")
        logger.info(F"Optimization method used is: {optimization_method}")

        if optimization_method.is_single_codon_optimization:
            result = optimize_sequence_by_single_codon(
                target_gene=module_input.sequence,
                organisms=module_input.organisms,
                optimization_cub_index=optimization_cub_index,
                optimization_method=optimization_method,
                tuning_param=module_input.tuning_parameter,
                skipped_codons_num=skipped_codons_num,
                run_summary=run_summary,
            )
            validate_module_output(original_sequence=module_input.sequence, new_sequence=result)
            return [result]

        if optimization_method.is_zscore_optimization:
            return ORFModule.optimize_sequence_by_zscore(
                module_input=module_input,
                optimization_cub_index=optimization_cub_index,
                optimization_method=optimization_method,
                skipped_codons_num=skipped_codons_num,
                run_summary=run_summary,
            )

        if optimization_method.is_single_organism_optimization:
            result = optimize_sequence_for_single_organism(
                target_gene=module_input.sequence,
                organisms=module_input.organisms,
                optimization_method=optimization_method,
                optimization_cub_index=optimization_cub_index,
                skipped_codons_num=skipped_codons_num,
                run_summary=run_summary,
            )
            validate_module_output(original_sequence=module_input.sequence, new_sequence=result)
            return [result]

        raise ValueError(F"optimization method {optimization_method} is invalid")

    @staticmethod
    def optimize_sequence_by_zscore(
            module_input: models.ModuleInput,
            optimization_cub_index: models.ORFOptimizationCubIndex,
            optimization_method: models.ORFOptimizationMethod,
            skipped_codons_num: int,
            run_summary: RunSummary,
    ):
        original_sequence = module_input.sequence
        target_genes = [original_sequence] + [
            synonymous_codon_permutation(seq=original_sequence, skipped_codons_num=skipped_codons_num) for _ in
            range(config["ORF"]["ZSCORE_INITIAL_PERMUTATIONS_NUM"])
        ]
        zscore_optimization = ORFModule.get_zscore_optimization_method(optimization_method)
        partial_zscore_optimization_method = partial(
            zscore_optimization,
            module_input=module_input,
            optimization_method=optimization_method,
            optimization_cub_index=optimization_cub_index,
            skipped_codons_num=skipped_codons_num,
            run_summary=run_summary,
        )

        results = []
        for target_gene in target_genes:
            logger.info(f"Running ORF optimization for initial sequence: {target_gene}")
            result = partial_zscore_optimization_method(sequence=target_gene)
            validate_module_output(original_sequence=original_sequence, new_sequence=result)
            results.append(result)

        return results

    @staticmethod
    def get_zscore_optimization_method(optimization_method: models.ORFOptimizationMethod) -> typing.Callable:
        if optimization_method.is_zscore_single_aa_optimization:
            return optimize_sequence_by_zscore_single_aa
        return optimize_sequence_by_zscore_bulk_aa
