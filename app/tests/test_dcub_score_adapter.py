import pytest

from modules import models
from modules.ORF.single_codon_optimization_method import _calculate_codons_loss
from modules.run_summary import RunSummary


def _organism(name, is_optimized, cai_weights):
    return models.Organism(
        name=name,
        is_optimized=is_optimized,
        optimization_priority=50,
        cai_profile=dict(cai_weights),
        tai_profile=dict(cai_weights),
        codon_frequencies={codon: 0.5 for codon in cai_weights},
    )


@pytest.fixture
def two_organisms():
    # Only the two-codon amino acids matter for these assertions; every codon
    # needs an entry so _get_max_organism_attribute_value never averages over
    # an empty profile.
    from modules.shared_functions_and_vars import nt_to_aa

    wanted_weights = {codon: 0.2 for codon in nt_to_aa}
    unwanted_weights = {codon: 0.2 for codon in nt_to_aa}
    # Make TGT clearly the wanted organism's favourite cysteine codon.
    wanted_weights["TGT"] = 0.9
    wanted_weights["TGC"] = 0.1
    unwanted_weights["TGT"] = 0.1
    unwanted_weights["TGC"] = 0.9
    return [
        _organism("wanted", True, wanted_weights),
        _organism("unwanted", False, unwanted_weights),
    ]


def test_calculate_codons_loss_can_run_twice_without_a_run_summary(two_organisms):
    """The adapter recomputes this table after ORFModule already wrote
    'orf_debug' into the run summary. Passing run_summary=None must skip the
    write rather than raising KeyError on the duplicate key."""
    run_summary = RunSummary()
    kwargs = dict(
        organisms=two_organisms,
        tuning_param=0.5,
        optimization_method=models.ORFOptimizationMethod.single_codon_diff,
        optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
    )

    first = _calculate_codons_loss(run_summary=run_summary, **kwargs)
    second = _calculate_codons_loss(run_summary=None, **kwargs)

    assert first == second
    assert "orf_debug" in run_summary.get()


def test_calculate_codons_loss_still_raises_on_a_genuine_duplicate_write(two_organisms):
    """Passing the same run_summary twice is still a bug and must still raise -
    run_summary=None is an opt-out, not a blanket suppression."""
    run_summary = RunSummary()
    kwargs = dict(
        organisms=two_organisms,
        tuning_param=0.5,
        optimization_method=models.ORFOptimizationMethod.single_codon_diff,
        optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
        run_summary=run_summary,
    )

    _calculate_codons_loss(**kwargs)
    with pytest.raises(KeyError):
        _calculate_codons_loss(**kwargs)


from modules.hotspot_avoidance.dcub_score_adapter import build_dcub_codon_table


def _module_input(organisms, optimization_method, sequence="ATG" * 20):
    return models.ModuleInput(
        organisms=organisms,
        sequence=sequence,
        output_path="",
        tuning_parameter=0.5,
        clusters_count=1,
        orf_optimization_method=optimization_method,
        orf_optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
    )


class TestSingleCodonFamily:
    def test_table_is_higher_is_better_and_argmax_matches_dcubs_own_choice(self, two_organisms):
        """DCUB's native table is a LOSS (lower is better) and its chosen codon
        is the argMIN. The adapter negates it, so the same codon must be the
        argMAX here - this is the property every downstream caller relies on."""
        method = models.ORFOptimizationMethod.single_codon_diff
        module_input = _module_input(two_organisms, method)

        loss_table = _calculate_codons_loss(
            organisms=two_organisms,
            tuning_param=0.5,
            optimization_method=method,
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            run_summary=None,
        )
        score_table = build_dcub_codon_table(
            module_input=module_input,
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence=module_input.sequence,
            skipped_codons_num=0,
        )

        for amino_acid, codon_losses in loss_table.items():
            best_by_loss = min(codon_losses, key=codon_losses.get)
            best_by_score = max(score_table[amino_acid], key=score_table[amino_acid].get)
            assert best_by_score == best_by_loss, f"disagreement for {amino_acid}"

    def test_every_codon_of_every_amino_acid_is_present(self, two_organisms):
        from modules.shared_functions_and_vars import synonymous_codons

        table = build_dcub_codon_table(
            module_input=_module_input(two_organisms, models.ORFOptimizationMethod.single_codon_diff),
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence="ATG" * 20,
            skipped_codons_num=0,
        )
        for amino_acid, codons in synonymous_codons.items():
            assert set(table[amino_acid]) == set(codons)


class TestSingleOrganismFamily:
    def test_argmax_matches_the_wanted_organisms_profile(self, two_organisms):
        from modules.ORF.single_organism_optimization_method import _get_optimal_codons

        method = models.ORFOptimizationMethod.single_wanted_organism
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index

        expected = _get_optimal_codons(
            organisms=two_organisms,
            optimization_method=method,
            optimization_cub_index=cub_index,
        )
        table = build_dcub_codon_table(
            module_input=_module_input(two_organisms, method),
            optimization_cub_index=cub_index,
            sequence="ATG" * 20,
            skipped_codons_num=0,
        )

        for amino_acid, expected_codon in expected.items():
            assert max(table[amino_acid], key=table[amino_acid].get) == expected_codon
