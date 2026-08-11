import tempfile
from pathlib import Path

import pytest

from modules import models
from modules.main import run_modules

EXAMPLE_DATA = Path(__file__).resolve().parent.parent.parent / "analysis" / "example_data"

# Captured by running this exact input through the current (pre-migration)
# pipeline. single_codon_diff is deterministic (no RNG involved), so this
# is exactly reproducible.
EXPECTED_FINAL_SEQUENCE = (
    "ATGGTTTCCAAGGGCGAGGAGGATAACATGGCTATCATTAAAGAGTTTATGCGTTTTAAAGTTCATATGGAAGGCTCTGTTAACGGCCATGAATTTGAA"
    "ATCGAAGGCGAAGGCGAAGGCCGTCCGTATGAAGGCACACAGACAGCTAAACTGAAAGTTACAAAAGGCGGCCCGCTGCCGTTTGCTTGGGATATCCTGT"
    "CTCCGCAGTTTATGTATGGCTCTAAAGCTTATGTTAAACATCCGGCTGATATCCCGGATTATCTGAAACTGTCTTTTCCGGAAGGCTTTAAATGGGAACG"
    "TGTTATGAACTTTGAAGATGGCGGCGTTGTTACAGTTACACAGGATTCTTCTCTGCAGGATGGCGAATTTATCTATAAAGTTAAACTGCGTGGCACAAAC"
    "TTTCCGTCTGATGGCCCGGTTATGCAGAAAAAAACAATGGGCTGGGAAGCTTCTTCTGAACGTATGTATCCGGAAGATGGCGCTCTGAAAGGCGAAATCA"
    "AACAGCGTCTGAAACTGAAAGATGGCGGCCATTATGATGCTGAAGTTAAAACAACATATAAAGCTAAAAAACCGGTTCAGCTGCCGGGCGCTTATAACGT"
    "TAACATCAAACTGGATATCACATCTCATAACGAAGATTATACAATCGTTGAACAGTATGAACGTGCTGAAGGCCGTCATTCTACAGGCGGCATGGATGA"
    "ACTGTATAAATAA"
)
EXPECTED_AVERAGE_DISTANCE_SCORE = 308.0068536646968
EXPECTED_WEAKEST_LINK_SCORE = 308.0068536646968
EXPECTED_RATIO_SCORE = 1.3201916787713533e-09


def _build_user_input(output_path: str) -> models.UserInput:
    mcherry_seq = (
        EXAMPLE_DATA.joinpath("mCherry_original.fasta")
        .read_text()
        .splitlines()[1]
        .strip()
    )
    return models.UserInput(
        sequence=mcherry_seq,
        tuning_param=50,
        organisms={
            "ecoli": models.OrganismRequest(
                genome_path=str(EXAMPLE_DATA / "Escherichia-coli.gb"),
                optimized=True,
                optimization_priority=50,
            ),
            "bsubtilis": models.OrganismRequest(
                genome_path=str(EXAMPLE_DATA / "Bacillus-subtilis.gb"),
                optimized=False,
                optimization_priority=50,
            ),
        },
        clusters_count=1,
        orf_optimization_method=models.ORFOptimizationMethod.single_codon_diff,
        orf_optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
        initiation_optimization_method=models.InitiationOptimizationMethod.original,
        output_path=output_path,
        evaluation_score="average_distance",
    )


def test_run_modules_end_to_end_mcherry_ecoli_vs_bsubtilis():
    with tempfile.TemporaryDirectory() as output_path:
        user_input = _build_user_input(output_path)
        result = run_modules(user_input, should_run_output_module=False)

    assert "error_message" not in result, result.get("error_message")
    final_evaluation = result["final_evaluation"]
    assert final_evaluation["final_sequence"] == EXPECTED_FINAL_SEQUENCE
    assert final_evaluation["average_distance_score"] == pytest.approx(EXPECTED_AVERAGE_DISTANCE_SCORE)
    assert final_evaluation["weakest_link_score"] == pytest.approx(EXPECTED_WEAKEST_LINK_SCORE)
    assert final_evaluation["ratio_score"] == pytest.approx(EXPECTED_RATIO_SCORE)
