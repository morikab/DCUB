import typing

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import validate_module_output

from .weak_folding_optimization import optimize_by_weak_folding

config = Configuration.get_config()
logger = LoggerFactory.get_logger()


class InitiationModule(object):
    @staticmethod
    def run_module(
            module_input: models.ModuleInput,
            run_summary: RunSummary,
    ) -> str:
        logger.info("\n##########################")
        logger.info("# Initiation #")
        logger.info("##########################")

        codons_num = config["INITIATION"]["NUMBER_OF_CODONS_TO_OPTIMIZE"]
        codons_num_in_nt = codons_num * 3
        logger.info(f"Running initiation optimization {module_input.initiation_optimization_method} on {codons_num} "
                    f"codons at the start of the ORF")
        if module_input.initiation_optimization_method == models.InitiationOptimizationMethod.original:
            logger.info(f"Taking {codons_num} codons from the start of the original ORF sequence:")
            optimized_sequence = module_input.sequence
            logger.info(f"Optimized sequence is: {optimized_sequence}")
        elif module_input.initiation_optimization_method == models.InitiationOptimizationMethod.external_module:
            logger.info(f"Taking optimized sequence from configuration value config.INITIATION.EXTERNAL_INITIATION_ORF")
            optimized_sequence = config["INITIATION"]["EXTERNAL_INITIATION_ORF"][:codons_num_in_nt]
            optimized_sequence += module_input.sequence[codons_num_in_nt:]
            logger.info(f"Optimized sequence is: {optimized_sequence}")
        elif module_input.initiation_optimization_method == models.InitiationOptimizationMethod.weak_folding:
            optimized_sequence = optimize_by_weak_folding(
                sequence=module_input.sequence,
                codons_num=codons_num,
                run_summary=run_summary,
            )
        else:
            raise ValueError(f"Initiation optimization method {module_input.initiation_optimization_method} is not "
                             f"supported")

        validate_module_output(original_sequence=module_input.sequence, new_sequence=optimized_sequence)
        return optimized_sequence
