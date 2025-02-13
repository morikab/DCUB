import typing
from collections import defaultdict

from numpy import average
from scipy.stats.mstats import gmean

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import change_all_codons_of_aa
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons
from modules.timer import Timer
from modules.ORF.calculating_cai import general_geomean


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
# In each round - check all single synonymous codon changes and calculate optimization score - take the best one
def optimize_sequence_by_zscore_single_aa(
        sequence: str,
        module_input: models.ModuleInput,
        optimization_cub_index: models.ORFOptimizationCubIndex,
        optimization_method: models.ORFOptimizationMethod,
        skipped_codons_num: int,
        run_summary: RunSummary,
        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"],
        should_dedup_codons: bool = config["ORF"]["DEDUP_CODONS"],
):
    """
    Iterative codon optimization:
    In each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter".
    """
    with Timer() as timer:
        initial_sequence = sequence
        previous_sequence_zscore = _calculate_zscore_for_sequence(
            sequence=sequence,
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            skipped_codons_num=skipped_codons_num,
        )
        initial_sequence_score = None
        previous_sequence_score = get_total_score(zscore=previous_sequence_zscore,
                                                  optimization_method=optimization_method,
                                                  tuning_parameter=module_input.tuning_parameter)

        aa_to_codon_mapping = defaultdict(str)
        iterations_count = 0
        iterations_summary = []
        for run in range(max_iterations):
            iterations_count = run + 1
            # Include also the sequence from the previous iteration
            sequence_to_zscore = {sequence: previous_sequence_zscore}
            tested_sequence_to_codon = defaultdict(list)
            for codon in nt_to_aa.keys():
                if nt_to_aa[codon] == "_" and optimization_cub_index.is_trna_adaptation_index:
                    # There is no point in optimizing stop codon by tAI weights, so keep the original codon
                    continue
                if len(synonymous_codons[nt_to_aa[codon]]) < 2:
                    # The codon does not have any other synonymous codons
                    continue
                tested_sequence = change_all_codons_of_aa(
                    seq=sequence,
                    selected_codon=codon,
                    skipped_codons_num=skipped_codons_num,
                )
                tested_sequence_to_codon[tested_sequence].append(codon)

                sequence_to_zscore[tested_sequence] = _calculate_zscore_for_sequence(
                    sequence=tested_sequence,
                    module_input=module_input,
                    optimization_cub_index=optimization_cub_index,
                    skipped_codons_num=skipped_codons_num,
                )

            if optimization_method.is_zscore_ratio_score_optimization:
                min_zscore = min(score.min_zscore for score in sequence_to_zscore.values())
                max_zscore = max(score.max_zscore for score in sequence_to_zscore.values())

                for score in sequence_to_zscore.values():
                    score.normalize(min_zscore=min_zscore, max_zscore=max_zscore)

            sequence_to_total_score = {
                sequence_option: get_total_score(zscore=sequence_score,
                                                 optimization_method=optimization_method,
                                                 tuning_parameter=module_input.tuning_parameter) for
                sequence_option, sequence_score in sequence_to_zscore.items()
            }
            if initial_sequence_score is None:
                initial_sequence_score = sequence_to_total_score[initial_sequence]

            new_sequence = max(sequence_to_total_score, key=sequence_to_total_score.get)
            selected_codons = tested_sequence_to_codon.get(new_sequence)
            best_new_sequence = new_sequence

            while len(sequence_to_total_score) > 0:
                non_valid_codon = None
                for codon in selected_codons:
                    frequencies = [o.codon_frequencies[codon] for o in module_input.organisms if o.is_optimized]
                    average_frequency = sum(frequencies) / len(frequencies)
                    if average_frequency < config["ORF"]["FREQUENCY_OPTIMIZATION_THRESHOLD"]:
                        logger.info(f"Skipping codon {codon} due to very low average frequency {average_frequency} in "
                                    f"wanted hosts.")
                        non_valid_codon = codon
                        break
                if non_valid_codon is None:
                    break
                sequence_to_total_score.pop(new_sequence)
                new_sequence = max(sequence_to_total_score, key=sequence_to_total_score.get)
                selected_codons = tested_sequence_to_codon.get(new_sequence)

            if len(sequence_to_total_score) < 1:
                new_sequence = best_new_sequence
                selected_codons = tested_sequence_to_codon.get(new_sequence)
                logger.info(f"Could not find codon that satisfies minimal average frequency in wanted "
                            f"hosts. Using the original optimal codons: {selected_codons}")

            new_sequence_score = sequence_to_total_score[new_sequence]

            # Summary information
            codon_to_score = {}
            for codon_sequence in sequence_to_total_score:
                for codon in tested_sequence_to_codon[codon_sequence]:
                    codon_to_score[codon] = sequence_to_total_score[codon_sequence]

            iteration_summary = {
                "selected_codons": [(selected_codon, nt_to_aa[selected_codon]) for selected_codon in selected_codons],
                "codon_to_score": codon_to_score,
                "sequence_score": new_sequence_score,
            }
            iterations_summary.append(iteration_summary)

            if new_sequence == sequence or previous_sequence_score > new_sequence_score:
                break

            sequence = new_sequence
            previous_sequence_score = new_sequence_score
            for selected_codon in selected_codons:
                # If the aa does not appear at all in the cds, this may give a faulty result (that should be consistent)
                aa_to_codon_mapping[nt_to_aa[selected_codon]] = selected_codon

    if should_dedup_codons:
        # Dedup repeating codons by using the second-best option of a synonymous codon for each selected codon
        for aa in aa_to_codon_mapping.keys():
            selected_aa_codon = None
            synonymous_codons_of_aa = synonymous_codons.get(aa).copy()
            for codon in synonymous_codons_of_aa:
                codon_index = sequence.find(codon)
                if codon_index != -1:
                    selected_aa_codon = codon
                    break
            synonymous_codons_of_aa.remove(selected_aa_codon)

            candidate_codons_sequences = []
            for candidate_codon in synonymous_codons_of_aa:
                # we determine the default codon based on the scores in the last iteration
                for seq in tested_sequence_to_codon.keys():
                    if candidate_codon in tested_sequence_to_codon[seq]:
                        candidate_codons_sequences.append(seq)
                        continue
            candidate_codons_scores = {seq: score for seq, score in
                                       sequence_to_total_score.items() if seq in candidate_codons_sequences}
            candidate_codon_entry = max(candidate_codons_scores, key=candidate_codons_scores.get)
            default_aa_codon = None
            for candidate_aa_default_codon in tested_sequence_to_codon.get(candidate_codon_entry):
                if nt_to_aa[candidate_aa_default_codon] == nt_to_aa[selected_aa_codon]:
                    default_aa_codon = candidate_aa_default_codon
                    break
            sequence = change_all_codons_of_aa(
                seq=sequence,
                selected_codon=selected_aa_codon,
                skipped_codons_num=skipped_codons_num,
                default_codon=default_aa_codon,
                should_dedup_codons=True,
            )

    orf_summary = {
        "initial_sequence": initial_sequence,
        "final_sequence": sequence,
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_codon_mapping,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": sequence_to_total_score[new_sequence],
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
    }
    run_summary.append_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def optimize_sequence_by_zscore_bulk_aa(
        sequence: str,
        module_input: models.ModuleInput,
        optimization_method: models.ORFOptimizationMethod,
        optimization_cub_index: models.ORFOptimizationCubIndex,
        skipped_codons_num: int,
        run_summary: RunSummary,
        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"],
        should_dedup_codons: bool = config["ORF"]["DEDUP_CODONS"],
):
    with Timer() as timer:
        initial_sequence = sequence
        initial_sequence_zscore = _calculate_zscore_for_sequence(
            sequence=sequence,
            module_input=module_input,
            optimization_cub_index=optimization_cub_index,
            skipped_codons_num=skipped_codons_num,
        )
        initial_sequence_score = None
        iterations_count = 0
        iterations_summary = []
        for run in range(max_iterations):
            logger.info(f"Starting run {run+1}.")
            iterations_count = run + 1
            codons_to_zscore = {}
            for codon in nt_to_aa.keys():
                candidate_codon_sequence = change_all_codons_of_aa(
                    seq=sequence,
                    selected_codon=codon,
                    skipped_codons_num=skipped_codons_num,
                )
                codons_to_zscore[codon] = _calculate_zscore_for_sequence(
                    sequence=candidate_codon_sequence,
                    module_input=module_input,
                    optimization_cub_index=optimization_cub_index,
                    skipped_codons_num=skipped_codons_num,
                )
            if optimization_method.is_zscore_ratio_score_optimization:
                min_zscore = min(score.min_zscore for score in codons_to_zscore.values())
                max_zscore = max(score.max_zscore for score in codons_to_zscore.values())

                for zscore in codons_to_zscore.values():
                    zscore.normalize(min_zscore=min_zscore, max_zscore=max_zscore)
                if initial_sequence_score is None:
                    initial_min_zscore = min(min_zscore, initial_sequence_zscore.min_zscore)
                    initial_max_zscore = max(max_zscore, initial_sequence_zscore.max_zscore)
                    initial_sequence_zscore.normalize(min_zscore=initial_min_zscore, max_zscore=initial_max_zscore)

            codons_to_total_score = {
                codon: get_total_score(zscore=zscore,
                                       optimization_method=optimization_method,
                                       tuning_parameter=module_input.tuning_parameter) for
                codon, zscore in codons_to_zscore.items()
            }

            if initial_sequence_score is None:
                initial_sequence_score = get_total_score(zscore=initial_sequence_zscore,
                                                         optimization_method=optimization_method,
                                                         tuning_parameter=module_input.tuning_parameter)
                previous_sequence_score = initial_sequence_score

            aa_to_selected_codon, aa_to_default_codon = _get_optimal_codon_per_aa(
                codons_to_total_score=codons_to_total_score,
                organisms=module_input.organisms,
            )

            # create new sequence by replacing all synonymous codons
            new_sequence = sequence
            for aa in aa_to_selected_codon:
                if aa == "_" and optimization_cub_index.is_trna_adaptation_index:
                    # There is no point in optimizing stop codon by tAI weights, so keep the original codon
                    continue
                new_sequence = change_all_codons_of_aa(
                    seq=new_sequence,
                    selected_codon=aa_to_selected_codon[aa],
                    default_codon=aa_to_default_codon[aa],
                    skipped_codons_num=skipped_codons_num,
                    should_dedup_codons=should_dedup_codons,
                )

            # Calculate score after all replacements
            zscore = _calculate_zscore_for_sequence(
                sequence=new_sequence,
                module_input=module_input,
                optimization_cub_index=optimization_cub_index,
                skipped_codons_num=skipped_codons_num,
            )

            if optimization_method.is_zscore_ratio_score_optimization:
                updated_min_zscore = min(min_zscore, zscore.min_zscore)
                updated_max_zscore = max(max_zscore, zscore.max_zscore)
                zscore.normalize(min_zscore=updated_min_zscore, max_zscore=updated_max_zscore)

            score = get_total_score(zscore=zscore,
                                    optimization_method=optimization_method,
                                    tuning_parameter=module_input.tuning_parameter)

            iteration_summary = {
                "aa_to_selected_codon": aa_to_selected_codon,
                "aa_to_default_codon": aa_to_default_codon,
                "sequence_score": score,
            }
            iterations_summary.append(iteration_summary)

            if new_sequence == sequence or previous_sequence_score > score:
                break
            else:
                sequence = new_sequence
                previous_sequence_score = score

    if new_sequence == sequence:
        aa_to_optimal_codon = aa_to_selected_codon
    elif len(iterations_summary) > 1:
        aa_to_optimal_codon = iterations_summary[-2]["aa_to_selected_codon"]
    else:
        aa_to_optimal_codon = []
        print(aa_to_optimal_codon)

    orf_summary = {
        "initial_sequence": initial_sequence,
        "final_sequence": sequence,
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_optimal_codon,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": previous_sequence_score,
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
    }

    run_summary.append_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def _get_optimal_codon_per_aa(
        codons_to_total_score: typing.Dict[str, float],
        organisms: typing.Sequence[models.Organism],
) -> typing.Tuple[typing.Dict[str, str], typing.Dict[str, str]]:
    aa_to_selected_codon = {}
    aa_to_default_codon = {}
    for aa in synonymous_codons.keys():
        codons_list = synonymous_codons[aa].copy()
        best_codon = _find_best_synonymous_codon_for_aa(codons_list=codons_list,
                                                        codons_to_score=codons_to_total_score)
        default_codon = None
        selected_codon = None
        while len(codons_list) > 0:
            candidate_codon = _find_best_synonymous_codon_for_aa(codons_list=codons_list,
                                                                 codons_to_score=codons_to_total_score)
            frequencies = [o.codon_frequencies[candidate_codon] for o in organisms if o.is_optimized]
            average_frequency = sum(frequencies) / len(frequencies)
            if average_frequency >= config["ORF"]["FREQUENCY_OPTIMIZATION_THRESHOLD"]:
                if selected_codon is None:
                    selected_codon = candidate_codon
                    logger.info(f"Found optimal codon {selected_codon} for amino acid {aa}. "
                                f"Now looking for second-best option.")
                    codons_list.remove(candidate_codon)
                else:
                    default_codon = candidate_codon
                    logger.info(f"Found default codon {default_codon} for amino acid {aa}.")
                    break
            else:
                logger.info(
                    f"Skipping codon {candidate_codon} due to very low average frequency {average_frequency} in "
                    f"wanted hosts.")
                codons_list.remove(candidate_codon)
        if selected_codon is None:
            logger.info(f"Could not find codon that satisfies minimal average frequency in wanted "
                        f"hosts. Using the codon with the best score: {best_codon}")
            selected_codon = best_codon

        aa_to_selected_codon[aa] = selected_codon
        aa_to_default_codon[aa] = default_codon
    return aa_to_selected_codon, aa_to_default_codon


# --------------------------------------------------------------
def _find_best_synonymous_codon_for_aa(codons_to_score: typing.Dict[str, float],
                                       codons_list: typing.Sequence[str]) -> str:
    aa_codons_to_score = {codon: score for codon, score in codons_to_score.items() if codon in codons_list}
    selected_aa_codon = max(aa_codons_to_score, key=aa_codons_to_score.get)
    return selected_aa_codon


# --------------------------------------------------------------
def _calculate_zscore_for_sequence(sequence: str,
                                   module_input: models.ModuleInput,
                                   optimization_cub_index: models.ORFOptimizationCubIndex,
                                   skipped_codons_num: int):
    sequence_for_calculation = sequence[skipped_codons_num*3:]
    optimization_cub_index_value = optimization_cub_index.value.lower()

    std_key = F"{optimization_cub_index_value}_std"
    average_key = F"{optimization_cub_index_value}_avg"
    weights = F"{optimization_cub_index_value}_profile"

    wanted_hosts_scores = []
    wanted_hosts_weights = []
    unwanted_hosts_scores = []
    unwanted_hosts_weights = []
    for organism in module_input.organisms:
        sigma = getattr(organism, std_key)
        miu = getattr(organism, average_key)
        profile = getattr(organism, weights)

        index_score = general_geomean([sequence_for_calculation], weights=profile)[0]
        organism_score = (index_score - miu) / sigma
        if organism.is_optimized:
            wanted_hosts_scores.append(organism_score)
            wanted_hosts_weights.append(organism.optimization_priority)
        else:
            unwanted_hosts_scores.append(organism_score)
            unwanted_hosts_weights.append(organism.optimization_priority)

    return models.SequenceZscores(
        initial_wanted_hosts_scores=wanted_hosts_scores,
        wanted_hosts_weights=wanted_hosts_weights,
        initial_unwanted_hosts_scores=unwanted_hosts_scores,
        unwanted_hosts_weights=unwanted_hosts_weights,
    )


# --------------------------------------------------------------
def get_total_score(zscore: models.SequenceZscores,
                    optimization_method: models.ORFOptimizationMethod,
                    tuning_parameter: float) -> float:
    if optimization_method.is_zscore_diff_score_optimization:
        return _calculate_zscore_diff_score(zscore=zscore,
                                            tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_ratio_score_optimization:
        return _calculate_zscore_ratio_score(zscore=zscore,
                                             tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_weakest_link_score_optimization:
        return _calculate_zscore_weakest_link_score(zscore=zscore,
                                                    tuning_parameter=tuning_parameter,
                                                    )


# --------------------------------------------------------------
def _calculate_zscore_diff_score(zscore: models.SequenceZscores,
                                 tuning_parameter: float) -> float:
    mean_opt_index = average(zscore.wanted_hosts_scores, weights=zscore.wanted_hosts_weights)
    mean_deopt_index = average(zscore.unwanted_hosts_scores, weights=zscore.unwanted_hosts_weights)
    return tuning_parameter * mean_opt_index - (1 - tuning_parameter) * mean_deopt_index


# --------------------------------------------------------------
def _calculate_zscore_ratio_score(zscore: models.SequenceZscores,
                                  tuning_parameter: float) -> float:
    mean_opt_index = gmean(zscore.wanted_hosts_scores, weights=zscore.wanted_hosts_weights)
    mean_deopt_index = gmean(zscore.unwanted_hosts_scores, weights=zscore.unwanted_hosts_weights)

    return (mean_opt_index ** tuning_parameter) / (mean_deopt_index ** (1 - tuning_parameter))


# --------------------------------------------------------------
def _calculate_zscore_weakest_link_score(zscore: models.SequenceZscores,
                                         tuning_parameter: float) -> float:
    weighted_optimized_organisms_scores = [zscore.wanted_hosts_scores[i] * zscore.wanted_hosts_weights[i] for i in
                                           range(len(zscore.wanted_hosts_scores))]
    weighted_deoptimized_organisms_scores = [
        zscore.unwanted_hosts_scores[i] * zscore.unwanted_hosts_weights[i] for
        i in range(len( zscore.unwanted_hosts_scores))
    ]

    return (tuning_parameter * min(weighted_optimized_organisms_scores) -
            (1 - tuning_parameter) * max(weighted_deoptimized_organisms_scores))
