import copy

import pytest

from modules import models
from modules.ORF.single_codon_optimization_method import _calculate_codons_loss
from modules.ORF.single_codon_optimization_method import orf_detailed_summary_key
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


def test_summary_key_is_scoped_by_cub_index():
    """A max_CAI_tAI run calls ORFModule twice against ONE RunSummary - once
    per index. An unscoped key made the second call collide, so max_CAI_tAI
    plus any single_codon_* method raised KeyError before this scoping."""
    cai = orf_detailed_summary_key(models.ORFOptimizationCubIndex.codon_adaptation_index)
    tai = orf_detailed_summary_key(models.ORFOptimizationCubIndex.trna_adaptation_index)

    assert cai == "orf_detailed_CAI"
    assert tai == "orf_detailed_tAI"
    assert cai != tai


def test_summary_key_is_scoped_by_stage():
    """Hotspot avoidance recomputes the same table. It records under its own
    stage key so a divergence between the two stages is visible rather than
    one silently overwriting the other."""
    index = models.ORFOptimizationCubIndex.codon_adaptation_index

    assert orf_detailed_summary_key(index) == "orf_detailed_CAI"
    assert (orf_detailed_summary_key(index, stage="hotspot_avoidance")
            == "hotspot_avoidance_detailed_CAI")


def test_both_cub_indexes_coexist_in_one_run_summary(two_organisms):
    """The regression this scoping exists for: two calls, one RunSummary,
    different indexes. Both entries must survive."""
    run_summary = RunSummary()
    kwargs = dict(
        organisms=two_organisms,
        tuning_param=0.5,
        optimization_method=models.ORFOptimizationMethod.single_codon_diff,
    )

    for index in (models.ORFOptimizationCubIndex.codon_adaptation_index,
                  models.ORFOptimizationCubIndex.trna_adaptation_index):
        _calculate_codons_loss(
            optimization_cub_index=index,
            run_summary=run_summary,
            summary_key=orf_detailed_summary_key(index),
            **kwargs,
        )

    summary = run_summary.get()
    assert "orf_detailed_CAI" in summary
    assert "orf_detailed_tAI" in summary
    # The two indexes read different profiles, so this also pins that the key
    # scoping preserved two DISTINCT tables rather than two copies of one.
    assert summary["orf_detailed_CAI"]["total_loss"]["C"] != {}
    assert set(summary["orf_detailed_CAI"]) == {
        "total_loss", "optimized_loss", "deoptimized_loss"
    }


def test_repeat_write_under_one_key_is_idempotent(two_organisms):
    """Hotspot avoidance patches each ORF candidate in turn, recomputing an
    identical table per candidate. Writing the same key twice must overwrite
    with the same value, not raise - the key is derived from exactly the
    inputs that determine the value."""
    run_summary = RunSummary()
    kwargs = dict(
        organisms=two_organisms,
        tuning_param=0.5,
        optimization_method=models.ORFOptimizationMethod.single_codon_diff,
        optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
        run_summary=run_summary,
        summary_key="orf_detailed_CAI",
    )

    first = _calculate_codons_loss(**kwargs)
    after_first = copy.deepcopy(run_summary.get()["orf_detailed_CAI"])
    second = _calculate_codons_loss(**kwargs)

    assert first == second
    assert run_summary.get()["orf_detailed_CAI"] == after_first


def test_max_cai_tai_single_codon_run_completes(two_organisms):
    """End-to-end regression for the key collisions. run_orf_optimization
    calls ORFModule twice for max_CAI_tAI - once per index - against one
    RunSummary. Both the per-codon loss dump and the "orf" summary used
    add_to_run_summary under an unscoped name, so the second call raised
    KeyError and this combination could not run at all."""
    from modules.main import run_orf_optimization

    module_input = models.ModuleInput(
        organisms=two_organisms,
        sequence="ATG" + "TGT" * 20,
        output_path="",
        tuning_parameter=0.5,
        clusters_count=1,
        orf_optimization_method=models.ORFOptimizationMethod.single_codon_diff,
        orf_optimization_cub_index=models.ORFOptimizationCubIndex.max_codon_trna_adaptation_index,
    )
    run_summary = RunSummary()

    cai, tai = run_orf_optimization(
        module_input=module_input, skipped_codons_num=0, run_summary=run_summary
    )

    assert len(cai) == 1 and len(tai) == 1
    summary = run_summary.get()
    assert {"orf_detailed_CAI", "orf_detailed_tAI"} <= set(summary)
    # "orf" appends, so both passes are retained rather than one overwriting
    # the other or raising.
    assert len(summary["orf"]) == 2


from modules.hotspot_avoidance.dcub_score_adapter import build_dcub_codon_table
from modules.hotspot_avoidance.dcub_score_adapter import build_dcub_score_fn
from modules.hotspot_avoidance.dcub_score_adapter import make_dcub_custom_score
from modules.hotspot_avoidance.dcub_score_adapter import rare_codons


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
            run_summary=RunSummary(),
            summary_key="orf_detailed_CAI",
        )
        score_table = build_dcub_codon_table(
            module_input=module_input,
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence=module_input.sequence,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        # _calculate_codons_loss returns each amino acid's codon dict sorted
        # ascending by loss, and _table_from_codon_loss's comprehension
        # preserves that key order while only negating the values. That
        # means min(codon_losses, ...) and max(score_table[aa], ...) would
        # both just return the dict's *first key* regardless of the actual
        # values - this test would pass even against an all-zeros score
        # table. Assert the values themselves first, so a degenerate table
        # cannot slip through silently:
        for amino_acid, codon_losses in loss_table.items():
            for codon, loss in codon_losses.items():
                assert score_table[amino_acid][codon] == pytest.approx(-loss)

        # Then restate the same claim in "argmax" vocabulary, the way the
        # docstring promises it. NOTE: this fixture ties the loss for most
        # amino acids - every codon but Cysteine's TGT/TGC is weighted
        # identically at 0.2 for both organisms - so min(codon_losses, ...)
        # is itself picking an arbitrary one of several equal-loss codons.
        # Re-sorting the score dict into alphabetical key order and taking a
        # second, independent max() (rather than checking whether the
        # loss table's own choice attains the maximum score) makes that
        # arbitrary tie-pick disagree purely on key-order grounds - verified
        # empirically: 18 of this fixture's 21 amino acids have a genuinely
        # tied minimum loss, and re-ordering flips the "winner" for all 18.
        # Checking that best_by_loss *attains* the max score - rather than
        # requiring it be THE unique argmax after an arbitrary reorder - is
        # the tie-safe version of the same check. It is also implied by the
        # per-value equality loop above (once every score is exactly -loss,
        # the argmin-of-loss trivially attains the max-of-score); it's kept
        # here for readability, matching the docstring's own vocabulary, not
        # as independent protection - the value-by-value loop above is what
        # defeats a degenerate/all-zeros table.
        for amino_acid, codon_losses in loss_table.items():
            best_by_loss = min(codon_losses, key=codon_losses.get)
            assert score_table[amino_acid][best_by_loss] == pytest.approx(
                max(score_table[amino_acid].values())
            ), f"disagreement for {amino_acid}"

    def test_every_codon_of_every_amino_acid_is_present(self, two_organisms):
        from modules.shared_functions_and_vars import synonymous_codons

        table = build_dcub_codon_table(
            module_input=_module_input(two_organisms, models.ORFOptimizationMethod.single_codon_diff),
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence="ATG" * 20,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )
        for amino_acid, codons in synonymous_codons.items():
            assert set(table[amino_acid]) == set(codons)


class TestSingleOrganismFamily:
    def test_argmax_matches_the_wanted_organisms_profile(self, two_organisms):
        from modules.ORF.single_organism_optimization_method import _get_optimal_codons

        # two_organisms gives every codon a uniform 0.2 weight except
        # TGT/TGC, so 20 of 21 amino acids here would only verify
        # tie-break-order-matching, not preference correctness. Rather than
        # modifying the shared fixture (Task 4's tests and Task 6 depend on
        # it as-is), build a local copy of the wanted organism's profile
        # with distinguishable per-codon weights across several amino
        # acids, so a real value difference - not ordering - decides the
        # answer for more than just Cysteine.
        wanted_organism, unwanted_organism = two_organisms
        modified_profile = dict(wanted_organism.cai_profile)
        modified_profile["GAT"] = 0.95
        modified_profile["GAC"] = 0.05
        modified_profile["AAG"] = 0.9
        modified_profile["AAA"] = 0.1
        modified_profile["TTT"] = 0.85
        modified_profile["TTC"] = 0.15
        organisms = [_organism("wanted", True, modified_profile), unwanted_organism]

        method = models.ORFOptimizationMethod.single_wanted_organism
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index

        expected = _get_optimal_codons(
            organisms=organisms,
            optimization_method=method,
            optimization_cub_index=cub_index,
        )
        table = build_dcub_codon_table(
            module_input=_module_input(organisms, method),
            optimization_cub_index=cub_index,
            sequence="ATG" * 20,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        for amino_acid, expected_codon in expected.items():
            assert max(table[amino_acid], key=table[amino_acid].get) == expected_codon


import math


def _zscore_organisms():
    """Organisms with the extra statistics the zscore path needs (cai_scores /
    tai_scores drive the cai_avg/cai_std properties used to standardize)."""
    from modules.shared_functions_and_vars import nt_to_aa

    wanted_weights = {codon: 0.2 for codon in nt_to_aa}
    wanted_weights["TGT"] = 0.9
    wanted_weights["TGC"] = 0.1
    unwanted_weights = {codon: 0.2 for codon in nt_to_aa}
    unwanted_weights["TGT"] = 0.1
    unwanted_weights["TGC"] = 0.9

    reference_scores = {f"gene_{index}": 0.1 * index for index in range(1, 11)}
    organisms = []
    for name, is_optimized, weights in (
        ("wanted", True, wanted_weights),
        ("unwanted", False, unwanted_weights),
    ):
        organism = models.Organism(
            name=name,
            is_optimized=is_optimized,
            optimization_priority=50,
            cai_profile=dict(weights),
            tai_profile=dict(weights),
            codon_frequencies={codon: 0.5 for codon in weights},
            cai_scores=dict(reference_scores),
            tai_scores=dict(reference_scores),
        )
        organisms.append(organism)
    return organisms


class TestZscoreFamilyIsScoredExactly:
    """The zscore family used to be approximated by a per-codon table built
    from "what if EVERY codon of this amino acid became X?". DCUB's real zscore
    optimization is positional and iterative, so that proxy measured a
    different quantity. It is now evaluated exactly, per trial sequence."""

    @pytest.mark.parametrize(
        "method",
        [
            models.ORFOptimizationMethod.zscore_bulk_aa_diff,
            models.ORFOptimizationMethod.zscore_bulk_aa_ratio,
            models.ORFOptimizationMethod.zscore_bulk_aa_weakest_link,
            models.ORFOptimizationMethod.zscore_single_aa_diff,
        ],
    )
    def test_score_is_finite(self, method):
        """Ratio methods geometric-mean the z-scores, which is undefined for the
        negative values standardization routinely produces - without fixed
        normalization bounds this silently returns nan."""
        organisms = _zscore_organisms()
        sequence = "ATGTGTTGCAAA" * 5
        score_fn = build_dcub_score_fn(
            module_input=_module_input(organisms, method, sequence=sequence),
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence=sequence,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        score = score_fn(sequence)
        assert isinstance(score, float)
        assert math.isfinite(score), f"{method} scored {score}"

    def test_score_matches_dcubs_own_evaluation_of_the_same_sequence(self):
        """The decisive property of exact scoring: for any candidate, the value
        must equal what DCUB's own zscore machinery computes for it - not an
        approximation that merely ranks similarly."""
        from modules.ORF.zscore_optimization_method import _calculate_zscore_for_sequence
        from modules.ORF.zscore_optimization_method import get_total_score

        method = models.ORFOptimizationMethod.zscore_bulk_aa_diff
        organisms = _zscore_organisms()
        sequence = "ATGTGTTGCAAA" * 5
        module_input = _module_input(organisms, method, sequence=sequence)
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index

        score_fn = build_dcub_score_fn(
            module_input=module_input,
            optimization_cub_index=cub_index,
            sequence=sequence,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        # A different candidate than the one the score fn was built around, so
        # this cannot pass by coincidence of being the reference sequence.
        candidate = sequence.replace("TGT", "TGC")
        expected = get_total_score(
            zscore=_calculate_zscore_for_sequence(
                sequence=candidate,
                module_input=module_input,
                optimization_cub_index=cub_index,
                skipped_codons_num=0,
            ),
            optimization_method=method,
            tuning_parameter=module_input.tuning_parameter,
        )
        assert score_fn(candidate) == pytest.approx(float(expected))

    def test_score_discriminates_between_candidates(self):
        """A scorer that returns the same number for every sequence would pass
        the finiteness check above while telling the optimizer nothing."""
        method = models.ORFOptimizationMethod.zscore_bulk_aa_diff
        organisms = _zscore_organisms()
        sequence = "ATGTGTTGCAAA" * 5
        score_fn = build_dcub_score_fn(
            module_input=_module_input(organisms, method, sequence=sequence),
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence=sequence,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        # TGT is the wanted host's favourite Cys codon and the unwanted host's
        # least favourite, so all-TGT must beat all-TGC for a diff method.
        assert score_fn(sequence.replace("TGC", "TGT")) > score_fn(sequence.replace("TGT", "TGC"))

    def test_ratio_normalization_bounds_are_fixed_across_calls(self):
        """Bounds are derived once and held. If they moved with each trial,
        successive scores would be incomparable - which is precisely what an
        optimizer must be able to do."""
        method = models.ORFOptimizationMethod.zscore_bulk_aa_ratio
        organisms = _zscore_organisms()
        sequence = "ATGTGTTGCAAA" * 5
        score_fn = build_dcub_score_fn(
            module_input=_module_input(organisms, method, sequence=sequence),
            optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            sequence=sequence,
            skipped_codons_num=0,
            run_summary=RunSummary(),
        )

        candidate = sequence.replace("TGT", "TGC")
        assert score_fn(candidate) == pytest.approx(score_fn(candidate))
        assert math.isfinite(score_fn(candidate))


def _frequency_organisms(wanted_frequencies):
    """Wanted/unwanted pair where TGT is the minimal-loss Cys codon but its
    frequency in the wanted host is caller-controlled."""
    from modules.shared_functions_and_vars import nt_to_aa

    wanted_weights = {codon: 0.2 for codon in nt_to_aa}
    unwanted_weights = {codon: 0.2 for codon in nt_to_aa}
    wanted_weights["TGT"], wanted_weights["TGC"] = 0.9, 0.1
    unwanted_weights["TGT"], unwanted_weights["TGC"] = 0.1, 0.9

    frequencies = {codon: 0.5 for codon in nt_to_aa}
    frequencies.update(wanted_frequencies)
    reference = {f"gene_{index}": 0.1 * index for index in range(1, 11)}
    return [
        models.Organism(name="wanted", is_optimized=True, optimization_priority=50,
                        cai_profile=dict(wanted_weights), tai_profile=dict(wanted_weights),
                        codon_frequencies=dict(frequencies),
                        cai_scores=dict(reference), tai_scores=dict(reference)),
        models.Organism(name="unwanted", is_optimized=False, optimization_priority=50,
                        cai_profile=dict(unwanted_weights), tai_profile=dict(unwanted_weights),
                        codon_frequencies={codon: 0.5 for codon in nt_to_aa},
                        cai_scores=dict(reference), tai_scores=dict(reference)),
    ]


class TestRareCodonFloor:
    """DCUB's _get_optimal_codon skips any codon whose mean frequency in the
    wanted hosts is below FREQUENCY_OPTIMIZATION_THRESHOLD. A plain
    negated-loss table knows nothing about that floor, so without this the
    optimizer could install a codon DCUB itself refused to select."""

    def test_rare_codons_are_identified(self):
        organisms = _frequency_organisms({"TGT": 0.02})
        assert "TGT" in rare_codons(organisms)
        assert "TGC" not in rare_codons(organisms)

    def test_amino_acids_with_no_compliant_codon_are_not_penalized(self):
        """_get_optimal_codon falls back to the minimal-loss codon when nothing
        clears the floor. Penalizing an unavoidable codon would steer nowhere
        and would skew scores between amino acids."""
        organisms = _frequency_organisms({"TGT": 0.02, "TGC": 0.01})
        penalized = rare_codons(organisms)
        assert "TGT" not in penalized
        assert "TGC" not in penalized

    def test_penalty_flips_the_argmax_to_dcubs_own_choice(self):
        """The regression this exists for, reproduced end to end: TGT has the
        minimal loss but is far too rare in the wanted host, so DCUB picks TGC.
        Without the penalty the score fn preferred TGT."""
        from modules.ORF.single_codon_optimization_method import _get_optimal_codon

        method = models.ORFOptimizationMethod.single_codon_diff
        organisms = _frequency_organisms({"TGT": 0.02})
        module_input = _module_input(organisms, method, sequence="ATG" + "TGT" * 10)
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index

        loss_table = _calculate_codons_loss(
            organisms=organisms, tuning_param=module_input.tuning_parameter,
            optimization_method=method, optimization_cub_index=cub_index,
            run_summary=RunSummary(), summary_key="probe",
        )
        dcub_choice = _get_optimal_codon(loss_table["C"].copy(), organisms)
        assert dcub_choice == "TGC", "fixture must make DCUB reject the minimal-loss codon"

        score_fn = build_dcub_score_fn(
            module_input=module_input, optimization_cub_index=cub_index,
            sequence=module_input.sequence, skipped_codons_num=0,
            run_summary=RunSummary(),
        )
        # Score one Cys codon in isolation, the choice the optimizer faces.
        assert score_fn("TGC") > score_fn("TGT")

        unpenalized = make_dcub_custom_score(build_dcub_codon_table(
            module_input=module_input, optimization_cub_index=cub_index,
            sequence=module_input.sequence, skipped_codons_num=0,
            run_summary=RunSummary(),
        ))
        assert unpenalized("TGT") > unpenalized("TGC"), (
            "without the floor the raw table prefers the rare codon - this is "
            "the divergence the penalty closes"
        )

    def test_floor_also_applies_to_the_zscore_family(self):
        """The floor is applied on top of the score, not by editing a table, so
        that it reaches the zscore family too - which has no table at all."""
        method = models.ORFOptimizationMethod.zscore_bulk_aa_diff
        organisms = _frequency_organisms({"TGT": 0.02})
        sequence = "ATGTGTTGCAAA" * 5
        module_input = _module_input(organisms, method, sequence=sequence)
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index

        score_fn = build_dcub_score_fn(
            module_input=module_input, optimization_cub_index=cub_index,
            sequence=sequence, skipped_codons_num=0, run_summary=RunSummary(),
        )
        all_tgt = sequence.replace("TGC", "TGT")
        all_tgc = sequence.replace("TGT", "TGC")

        # Exactly the opposite of test_score_discriminates_between_candidates,
        # which uses the same weights but compliant frequencies: there all-TGT
        # wins on the zscore alone; here the floor overrides it.
        assert score_fn(all_tgc) > score_fn(all_tgt)

    def test_no_penalty_when_every_codon_clears_the_floor(self):
        organisms = _frequency_organisms({})
        assert rare_codons(organisms) == frozenset()

        method = models.ORFOptimizationMethod.single_codon_diff
        module_input = _module_input(organisms, method, sequence="ATG" + "TGT" * 10)
        cub_index = models.ORFOptimizationCubIndex.codon_adaptation_index
        table = build_dcub_codon_table(
            module_input=module_input, optimization_cub_index=cub_index,
            sequence=module_input.sequence, skipped_codons_num=0,
            run_summary=RunSummary(),
        )
        score_fn = build_dcub_score_fn(
            module_input=module_input, optimization_cub_index=cub_index,
            sequence=module_input.sequence, skipped_codons_num=0,
            run_summary=RunSummary(),
        )
        assert score_fn("TGT") == pytest.approx(make_dcub_custom_score(table)("TGT"))


def test_get_optimal_codon_falls_back_instead_of_crashing():
    """`.keys()[0]` raises TypeError on Python 3, so every amino acid whose
    synonymous codons are ALL below the floor crashed the run rather than
    falling back to the minimal-loss codon as the log message promises."""
    from modules.ORF.single_codon_optimization_method import _get_optimal_codon

    organisms = _frequency_organisms({"TGT": 0.02, "TGC": 0.01})
    # Ordered ascending by loss, the order _calculate_codons_loss emits.
    candidates = {"TGT": 0.005, "TGC": 0.725}
    assert _get_optimal_codon(candidates, organisms) == "TGT"
