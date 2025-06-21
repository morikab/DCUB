import typing

import numpy as np
from numpy import average
from scipy.stats.mstats import gmean

from logger_factory.logger_factory import LoggerFactory
from modules import models as main_models
from modules import shared_functions_and_vars
from modules.run_summary import RunSummary
from modules.ORF.calculating_cai import general_geomean

from . import models

logger = LoggerFactory.get_logger()


class EvaluationModule(object):
    @staticmethod
    def run_module(final_sequence: str,
                   module_input: main_models.ModuleInput,
                   optimization_cub_index: main_models.ORFOptimizationCubIndex,
                   skipped_codons_num: int,
                   run_summary: RunSummary) -> models.EvaluationModuleResult:
        optimization_cub_index_value = optimization_cub_index.value.lower()
        initial_sequence = module_input.sequence
        std = f"{optimization_cub_index_value}_std"
        weights = f"{optimization_cub_index_value}_profile"

        optimized_organisms_zscores = []
        optimized_organisms_weights = []
        deoptimized_organisms_zscores = []
        deoptimized_organisms_weights = []
        zscores_for_normalization = []

        # Consider only relevant part of the initial and final sequences
        initial_sequence_for_evaluation = initial_sequence[skipped_codons_num * 3:]
        final_sequence_for_evaluation = final_sequence[skipped_codons_num * 3:]
        logger.info(f"Calculating evaluation score based on initial sequence: {initial_sequence_for_evaluation} and "
                    f"final sequence: {final_sequence_for_evaluation}")

        organisms_evaluation_summary = []
        for organism in module_input.organisms:
            sigma = getattr(organism, std)
            profile = getattr(organism, weights)
            edge_case_sequences = EvaluationModule._get_sequences_for_normalization(
                sequence=initial_sequence_for_evaluation,
                weights=profile,
            )

            cub_scores = general_geomean(
                [initial_sequence_for_evaluation, final_sequence_for_evaluation, *edge_case_sequences],
                weights=profile,
            )
            initial_score = cub_scores[0]
            final_score = cub_scores[1]
            edge_case_scores = cub_scores[2:]
            organism_zscore = (final_score - initial_score) / sigma
            zscores_for_normalization.extend([(score - initial_score) / sigma for score in edge_case_scores])
            if organism.is_optimized:
                optimized_organisms_zscores.append(organism_zscore)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_zscores.append(organism_zscore)
                deoptimized_organisms_weights.append(organism.optimization_priority)

            organism_summary = {
                "name": organism.name,
                "is_wanted": organism.is_optimized,
                F"{optimization_cub_index_value}_initial_score": initial_score,
                F"{optimization_cub_index_value}_final_score": final_score,
                "dist_score": organism_zscore,
            }
            organisms_evaluation_summary.append(organism_summary)

        average_distance_score = EvaluationModule._calculate_average_distance_score(
            optimized_organisms_scores=optimized_organisms_zscores,
            deoptimized_organisms_scores=deoptimized_organisms_zscores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=module_input.tuning_parameter,
        )

        weakest_link_score = EvaluationModule._calculate_weakest_link_score(
            optimized_organisms_scores=optimized_organisms_zscores,
            deoptimized_organisms_scores=deoptimized_organisms_zscores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=module_input.tuning_parameter,
        )

        ratio_score = EvaluationModule._calculate_ratio_score(
            optimized_organisms_scores=optimized_organisms_zscores,
            deoptimized_organisms_scores=deoptimized_organisms_zscores,
            scores_for_normalization=zscores_for_normalization,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=module_input.tuning_parameter,
        )

        evaluation_result = models.EvaluationModuleResult(
            sequence=final_sequence,
            average_distance_score=average_distance_score,
            weakest_link_score=weakest_link_score,
            ratio_score=ratio_score,
        )

        evaluation_summary = {
            "organisms": organisms_evaluation_summary,
            **evaluation_result.summary,
            "init_evaluation_sequence": initial_sequence_for_evaluation,
            "final_evaluation_sequence": final_sequence_for_evaluation,
        }
        run_summary.append_to_run_summary("evaluation", evaluation_summary)

        return evaluation_result

    # --------------------------------------------------------------
    @staticmethod
    def _get_sequences_for_normalization(sequence: str, weights: typing.Dict[str, float]) -> typing.Sequence[str]:
        def _get_sequence_per_aggregation_method(aggregation_method):
            aa_to_codon = {}
            for aa, aa_codons in shared_functions_and_vars.synonymous_codons.items():
                codon_weights = {codon: weights[codon] for codon in aa_codons if weights.get(codon)}
                if not codon_weights:
                    continue
                aa_to_codon[aa] = aggregation_method(codon_weights, key=codon_weights.get)
            final_sequence = sequence
            for aa in aa_to_codon.keys():
                final_sequence = shared_functions_and_vars.change_all_codons_of_aa(
                    seq=final_sequence,
                    selected_codon=aa_to_codon[aa],
                )
            return final_sequence

        return [_get_sequence_per_aggregation_method(method) for method in (min, max)]

    # --------------------------------------------------------------
    @staticmethod
    def _calculate_average_distance_score(optimized_organisms_scores: typing.Sequence[float],
                                          deoptimized_organisms_scores: typing.Sequence[float],
                                          optimized_organisms_weights: typing.Sequence[float],
                                          deoptimized_organisms_weights: typing.Sequence[float],
                                          tuning_parameter: float) -> float:
        mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
        mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
        return tuning_parameter * mean_opt_index - (1-tuning_parameter) * mean_deopt_index

    # --------------------------------------------------------------
    @staticmethod
    def _calculate_ratio_score(optimized_organisms_scores: typing.Sequence[float],
                               deoptimized_organisms_scores: typing.Sequence[float],
                               scores_for_normalization: typing.Sequence[float],
                               optimized_organisms_weights: typing.Sequence[float],
                               deoptimized_organisms_weights: typing.Sequence[float],
                               tuning_parameter: float) -> float:
        all_scores = np.array([*scores_for_normalization, *optimized_organisms_scores, *deoptimized_organisms_scores])
        min_score = all_scores.min()
        max_score = all_scores.max()
        epsilon = 1e-6
        normalized_optimized_organisms_scores = [
            (score - min_score + epsilon) / (max_score - min_score + epsilon) for
            score in optimized_organisms_scores
        ]
        normalized_deoptimized_organisms_scores = [
            (score - min_score + epsilon) / (max_score - min_score + epsilon) for
            score in deoptimized_organisms_scores
        ]
        mean_opt_index = gmean(normalized_optimized_organisms_scores, weights=optimized_organisms_weights)
        mean_deopt_index = gmean(normalized_deoptimized_organisms_scores, weights=deoptimized_organisms_weights)

        return (mean_opt_index ** tuning_parameter) / (mean_deopt_index ** (1 - tuning_parameter))
    # --------------------------------------------------------------
    @staticmethod
    def _calculate_weakest_link_score(optimized_organisms_scores: typing.Sequence[float],
                                      deoptimized_organisms_scores: typing.Sequence[float],
                                      optimized_organisms_weights: typing.Sequence[float],
                                      deoptimized_organisms_weights: typing.Sequence[float],
                                      tuning_parameter: float) -> float:
        weighted_optimized_organisms_scores = [
            optimized_organisms_scores[i] * optimized_organisms_weights[i] for i in
            range(len(optimized_organisms_scores))
        ]
        weighted_deoptimized_organisms_scores = [
            deoptimized_organisms_scores[i] * deoptimized_organisms_weights[i] for i in
            range(len(deoptimized_organisms_scores))
        ]
        return tuning_parameter * min(weighted_optimized_organisms_scores) - (1-tuning_parameter) * max(
            weighted_deoptimized_organisms_scores
        )
