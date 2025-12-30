import os
import shutil
import time
import traceback
from pathlib import Path
import typing

from logger_factory.logger_factory import LoggerFactory
from modules.configuration import Configuration
from modules.run_summary import RunSummary

# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
# if artifacts_directory.exists() and artifacts_directory.is_dir():
#     shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)
artifacts_path = str(artifacts_directory)

from modules import initiation
from modules import ORF
from modules import sequence_family
from modules import user_IO
from modules.evaluation import models as evaluation_models
from modules.evaluation.evaluation import EvaluationModule
from modules import models

logger = LoggerFactory.get_logger()
config = Configuration.read_config()


def run_input_processing(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None) -> models.ModuleInput:
    run_summary = RunSummary()
    return user_IO.UserInputModule.run_module(user_inp_raw=user_input_dict, run_summary=run_summary)


def run_modules(user_input_dict: typing.Dict[str, typing.Any],
                should_run_output_module: bool = True) -> typing.Dict[str, typing.Any]:
    run_summary = RunSummary()
    final_output = {}
    try:
        before_parsing_input = time.time()
        initiation_optimized_codons_num = config["INITIATION"]["NUMBER_OF_CODONS_TO_OPTIMIZE"]
        module_input = user_IO.UserInputModule.run_module(
            user_inp_raw=user_input_dict,
            initiation_optimized_codons_num=initiation_optimized_codons_num,
            run_summary=run_summary,
        )
        after_parsing_input = time.time()
        logger.info(F"Total input processing time: {after_parsing_input-before_parsing_input}")
        # ####################################### Initiation Optimization #################################
        initiation_optimized_sequence = initiation.InitiationModule.run_module(
            module_input=module_input,
            run_summary=run_summary,
        )
        # ####################################### ORF Optimization ########################################
        module_input.sequence = initiation_optimized_sequence
        clustered_module_inputs = sequence_family.SequenceFamilyModule.run_module(module_input)
        clustered_orf_results = []
        for module_input_cluster in clustered_module_inputs:
            clustered_orf_results.append(run_orf_optimization(
                module_input=module_input_cluster,
                skipped_codons_num=initiation_optimized_codons_num,
                run_summary=run_summary,
            ))
        # ####################################### Evaluation ##############################################
        evaluation_results = []
        for cds_nt_final_cai, cds_nt_final_tai in clustered_orf_results:
            evaluation_result = run_evaluation(
                module_input=module_input,
                cds_nt_final_cai=cds_nt_final_cai,
                cds_nt_final_tai=cds_nt_final_tai,
                skipped_codons_num=initiation_optimized_codons_num,
                run_summary=run_summary,
            )
            evaluation_results.append(evaluation_result)
        # ###################################### Output Handling ##########################################
        output_path = module_input.output_path or str(artifacts_directory)
        run_summary.save_run_summary(output_path)

        final_output = run_summary.get()

        # TODO - handle multiple results in output generation module
        evaluation_result = evaluation_results[0]
        if should_run_output_module:
            final_output["zip_output_file_path"] = user_IO.UserOutputModule.run_module(
                cds_sequence=evaluation_result.sequence,
                output_path=output_path,
                artifacts_path=artifacts_path,
            )
    except:
        logger.error("Encountered unknown error when running modules.")
        exception_str = traceback.format_exc()
        final_output = {
            "error_message": exception_str,
        }
    finally:
        logger.info("Final output: %s", final_output)

    return final_output


def run_evaluation(
    module_input: models.ModuleInput,
    cds_nt_final_cai: typing.Sequence[str],
    cds_nt_final_tai: typing.Sequence[str],
    skipped_codons_num: int,
    run_summary: RunSummary,
) -> evaluation_models.EvaluationModuleResult:
    logger.info('\n##########################')
    logger.info('# EVALUATION #')
    logger.info('##########################')
    tai_evaluation_results = [EvaluationModule.run_module(
        final_sequence=cds_nt_tai,
        module_input=module_input,
        optimization_cub_index=models.ORFOptimizationCubIndex.trna_adaptation_index,
        skipped_codons_num=skipped_codons_num,
        run_summary=run_summary,
    ) for cds_nt_tai in cds_nt_final_tai]

    cai_evaluation_results = [EvaluationModule.run_module(
        final_sequence=cds_nt_cai,
        module_input=module_input,
        optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
        skipped_codons_num=skipped_codons_num,
        run_summary=run_summary,
    ) for cds_nt_cai in cds_nt_final_cai]

    evaluation_result = choose_orf_optimization_result(
        tai_evaluation_results=tai_evaluation_results,
        cai_evaluation_results=cai_evaluation_results,
        evaluation_score=module_input.evaluation_score,
    )
    logger.info(f"Final evaluation result: {evaluation_result.summary}")
    run_summary.add_to_run_summary("final_evaluation", evaluation_result.summary)
    return evaluation_result


def choose_orf_optimization_result(
        evaluation_score: models.EvaluationScore,
        tai_evaluation_results: typing.Optional[typing.Sequence[evaluation_models.EvaluationModuleResult]],
        cai_evaluation_results: typing.Optional[typing.Sequence[evaluation_models.EvaluationModuleResult]],
) -> evaluation_models.EvaluationModuleResult:
    tai_evaluation_results = tai_evaluation_results or []
    cai_evaluation_results = cai_evaluation_results or []
    evaluation_results = tai_evaluation_results + cai_evaluation_results

    best_evaluation_result = None
    for evaluation_result in evaluation_results:
        logger.info(F"Analyzing evaluation result: {evaluation_result.summary}")
        if best_evaluation_result is None:
            best_evaluation_result = evaluation_result
            continue
        if evaluation_result.get_score(evaluation_score) > best_evaluation_result.get_score(evaluation_score):
            best_evaluation_result = evaluation_result

    return best_evaluation_result


def run_orf_optimization(
        module_input: models.ModuleInput,
        skipped_codons_num: int,
        run_summary: RunSummary,
) -> typing.Tuple[typing.Sequence[str], typing.Sequence[str]]:
    optimization_cub_index = module_input.orf_optimization_cub_index
    optimization_method = module_input.orf_optimization_method
    cds_nt_final_cai = []
    cds_nt_final_tai = []
    if optimization_cub_index.is_trna_adaptation_index:
        cds_nt_final_tai = ORF.ORFModule.run_module(
            module_input=module_input,
            optimization_cub_index=models.ORFOptimizationCubIndex.trna_adaptation_index,
            optimization_method=optimization_method,
            skipped_codons_num=skipped_codons_num,
            run_summary=run_summary,
        )

    if optimization_cub_index.is_codon_adaptation_index:
        cds_nt_final_cai = ORF.ORFModule.run_module(
            module_input=module_input,
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            optimization_method=optimization_method,
            skipped_codons_num=skipped_codons_num,
            run_summary=run_summary,
        )

    return cds_nt_final_cai, cds_nt_final_tai


def run_orf_module(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]]):
    run_summary = RunSummary()
    module_input = user_IO.UserInputModule.run_module(user_inp_raw=user_input_dict, run_summary=run_summary)

    orf_optimization_cub_index = module_input.orf_optimization_cub_index
    optimization_method = module_input.orf_optimization_method

    skipped_codons_num = config["INITIATION"]["NUMBER_OF_CODONS_TO_OPTIMIZE"]
    if orf_optimization_cub_index.is_trna_adaptation_index:
        logger.info("tAI information:")
        trna_adaptation_index = orf_optimization_cub_index.trna_adaptation_index
        return ORF.ORFModule.run_module(module_input=module_input,
                                        optimization_cub_index=trna_adaptation_index,
                                        optimization_method=optimization_method,
                                        skipped_codons_num=skipped_codons_num,
                                        run_summary=run_summary)

    codon_adaptation_index = orf_optimization_cub_index.codon_adaptation_index
    return ORF.ORFModule.run_module(module_input=module_input,
                                    optimization_cub_index=codon_adaptation_index,
                                    optimization_method=optimization_method,
                                    skipped_codons_num=skipped_codons_num,
                                    run_summary=run_summary)
