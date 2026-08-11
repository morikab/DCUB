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
