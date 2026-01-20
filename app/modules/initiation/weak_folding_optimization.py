# FIXME
# from RNA import RNA

from logger_factory.logger_factory import LoggerFactory
from modules import shared_functions_and_vars
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.timer import Timer

logger = LoggerFactory.get_logger()
config = Configuration.get_config()


def _calculate_folding_strength_for_sequence(
    sequence: str,
):
    """
    MFE (minimal free energy) of the sequence is returned in units of kcal/mol
    """
    from RNA import RNA
    fold_compound = RNA.fold_compound(sequence)
    _, mfe = fold_compound.mfe()
    return mfe


# TODO - need to fix with Shimshi's code
def optimize_by_weak_folding(
        sequence: str,
        codons_num: int,
        run_summary: RunSummary,
        max_iterations: int = config["INITIATION"]["MAX_ITERATIONS"],
) -> str:
    codon_size = 3
    prefix_size_in_nt = codons_num * codon_size
    with Timer() as timer:
        initial_prefix = sequence[:prefix_size_in_nt]
        logger.info(f"Initial initiation sequence: \n {initial_prefix}")
        prefix = initial_prefix
        previous_mfe = _calculate_folding_strength_for_sequence(
            sequence=prefix,
        )
        initial_mfe = None
        iterations_count = 0
        # Single codon replacement
        for run in range(max_iterations):
            iterations_count = run + 1
            # Include also the sequence from the previous iteration
            sequence_to_mfe = {prefix: previous_mfe}
            for i in range(0, len(prefix), codon_size):
                codon = prefix[i:i+3]
                synonymous_codons = shared_functions_and_vars.synonymous_codons[
                    shared_functions_and_vars.nt_to_aa[codon]
                ]
                for synonymous_codon in synonymous_codons:
                    if synonymous_codon != codon:
                        candidate_prefix = prefix[:i] + synonymous_codon + prefix[i+codon_size:]
                        sequence_to_mfe[candidate_prefix] = _calculate_folding_strength_for_sequence(
                            sequence=candidate_prefix,
                        )

            if initial_mfe is None:
                initial_mfe = sequence_to_mfe[initial_prefix]

            # The more negative is the mfe - the stronger the folding energy of the mRNA.
            # We look for maximal mfe value as we're optimizing for a WEAK folding energy.
            new_prefix = max(sequence_to_mfe, key=sequence_to_mfe.get)

            if new_prefix == prefix:
                break
            else:
                prefix = new_prefix
                previous_mfe = sequence_to_mfe[prefix]

    logger.info(f"Final initiation sequence: \n {prefix}")
    new_sequence = prefix + sequence[prefix_size_in_nt:]
    initiation_summary = {
        "initial_sequence": sequence,
        "final_sequence": new_sequence,
        "prefix_codons_num": codons_num,
        "iterations_count": iterations_count,
        "initial_mfe": initial_mfe,
        "final_mfe": sequence_to_mfe[prefix],
        "run_time": timer.elapsed_time,
    }
    run_summary.append_to_run_summary("initiation", initiation_summary)

    return new_sequence
