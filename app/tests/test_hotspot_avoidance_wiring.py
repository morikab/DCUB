from modules import models


class TestEnableHotspotAvoidanceFlag:
    def test_module_input_defaults_to_disabled(self):
        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG",
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
        )
        assert module_input.enable_hotspot_avoidance is False

    def test_user_input_defaults_to_disabled(self):
        user_input = models.UserInput(
            sequence="ATG",
            tuning_param=0.5,
            organisms={},
            clusters_count=1,
            orf_optimization_method=models.ORFOptimizationMethod.single_codon_diff,
            orf_optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            initiation_optimization_method=models.InitiationOptimizationMethod.original,
            output_path="",
            evaluation_score="average_distance",
        )
        assert user_input.enable_hotspot_avoidance is False

    def test_user_input_accepts_the_flag(self):
        user_input = models.UserInput(
            sequence="ATG",
            tuning_param=0.5,
            organisms={},
            clusters_count=1,
            orf_optimization_method=models.ORFOptimizationMethod.single_codon_diff,
            orf_optimization_cub_index=models.ORFOptimizationCubIndex.codon_adaptation_index,
            initiation_optimization_method=models.InitiationOptimizationMethod.original,
            output_path="",
            evaluation_score="average_distance",
            enable_hotspot_avoidance=True,
        )
        assert user_input.enable_hotspot_avoidance is True

    def test_module_input_summary_reports_the_flag(self):
        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG",
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
            enable_hotspot_avoidance=True,
        )
        assert module_input.summary["enable_hotspot_avoidance"] is True


class TestDedupRetirement:
    def test_dedup_is_forced_off_when_hotspot_avoidance_is_enabled(self):
        """DCUB's repeat-avoidance heuristic and ESO's slippage/recombination
        detection solve the same problem from different angles. When ESO is
        doing it properly, DCUB's heuristic is retired for that run."""
        from modules.ORF.orf_main import ORFModule

        assert ORFModule.should_dedup_codons(enable_hotspot_avoidance=True) is False

    def test_dedup_follows_config_when_hotspot_avoidance_is_disabled(self, monkeypatch):
        # config["ORF"]["DEDUP_CODONS"] is False today, so asserting against
        # it directly cannot tell the real helper apart from one that
        # ignores the flag and always returns False. Drive the True branch
        # explicitly so the "retire dedup" behaviour is actually exercised.
        from modules.configuration import Configuration
        from modules.ORF import orf_main

        expected = Configuration.get_config()["ORF"]["DEDUP_CODONS"]
        assert orf_main.ORFModule.should_dedup_codons(enable_hotspot_avoidance=False) == expected

        monkeypatch.setitem(orf_main.config["ORF"], "DEDUP_CODONS", True)
        assert orf_main.ORFModule.should_dedup_codons(enable_hotspot_avoidance=False) is True
        assert orf_main.ORFModule.should_dedup_codons(enable_hotspot_avoidance=True) is False


class TestRunHotspotAvoidance:
    def test_disabled_returns_candidates_untouched(self):
        from modules.main import run_hotspot_avoidance

        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG" * 10,
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
            enable_hotspot_avoidance=False,
        )
        cai = ["AAA" * 10]
        tai = ["CCC" * 10]

        patched_cai, patched_tai, summaries = run_hotspot_avoidance(
            module_input=module_input,
            cds_nt_final_cai=cai,
            cds_nt_final_tai=tai,
            skipped_codons_num=0,
        )

        assert patched_cai == cai
        assert patched_tai == tai
        assert summaries == {}

    def test_enabled_patches_every_candidate_in_both_lists(self, monkeypatch):
        """Every candidate must be patched before evaluation, not just one -
        selection IS the evaluation step."""
        from modules import main as main_module
        from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult

        seen = []

        def fake_run_module(sequence, module_input, optimization_cub_index, skipped_codons_num):
            seen.append((sequence, optimization_cub_index))
            return HotspotPatchResult(
                sequence_before=sequence,
                sequence_after=sequence.replace("AAA", "AAG", 1),
                num_edits=1,
                detected_sites={"recombination": 0, "slippage": 1, "motifs": 0},
                warnings=[],
            )

        monkeypatch.setattr(
            main_module.HotspotAvoidanceModule, "run_module", staticmethod(fake_run_module)
        )

        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG" * 10,
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
            enable_hotspot_avoidance=True,
        )
        cai = ["AAA" + "CCC" * 9, "AAA" + "GGG" * 9]
        tai = ["AAA" + "TTT" * 9]

        patched_cai, patched_tai, summaries = main_module.run_hotspot_avoidance(
            module_input=module_input,
            cds_nt_final_cai=cai,
            cds_nt_final_tai=tai,
            skipped_codons_num=0,
        )

        assert [sequence for sequence, _ in seen] == cai + tai, (
            "every candidate in both lists must be patched"
        )
        # The CAI list must be scored with codon_adaptation_index and the tAI list with
        # trna_adaptation_index - never max_codon_trna_adaptation_index. That value would
        # resolve to a "max_cai_tai_profile" attribute Organism does not have, so
        # build_dcub_codon_table's getattr would fall back to {} and hand the optimizer an
        # all-zero codon table with no error raised anywhere.
        assert [index for _, index in seen] == (
            [models.ORFOptimizationCubIndex.codon_adaptation_index] * len(cai)
            + [models.ORFOptimizationCubIndex.trna_adaptation_index] * len(tai)
        )
        assert patched_cai == [candidate.replace("AAA", "AAG", 1) for candidate in cai]
        assert patched_tai == [candidate.replace("AAA", "AAG", 1) for candidate in tai]
        # The summary lookup is keyed by the PATCHED sequence, so the winner can
        # be resolved after evaluation picks it.
        for patched in patched_cai + patched_tai:
            assert summaries[patched]["num_edits"] == 1

    def test_translation_preserving_failure_is_caught_by_validate_module_output(self, monkeypatch):
        from modules import main as main_module
        from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult

        def broken_run_module(sequence, module_input, optimization_cub_index, skipped_codons_num):
            return HotspotPatchResult(
                sequence_before=sequence,
                sequence_after="ATG",  # wrong length
                num_edits=1,
                detected_sites={"recombination": 0, "slippage": 0, "motifs": 0},
                warnings=[],
            )

        monkeypatch.setattr(
            main_module.HotspotAvoidanceModule, "run_module", staticmethod(broken_run_module)
        )

        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG" * 10,
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
            enable_hotspot_avoidance=True,
        )

        import pytest

        with pytest.raises(RuntimeError):
            main_module.run_hotspot_avoidance(
                module_input=module_input,
                cds_nt_final_cai=["AAA" * 10],
                cds_nt_final_tai=[],
                skipped_codons_num=0,
            )
