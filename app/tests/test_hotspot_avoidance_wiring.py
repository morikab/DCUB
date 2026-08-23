from modules import models
from modules.run_summary import RunSummary


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
            run_summary=RunSummary(),
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

        def fake_run_module(sequence, module_input, optimization_cub_index, skipped_codons_num,
                            run_summary, compute_motifs=None):
            # run_summary is captured, not ignored: the codon-loss table the
            # real module builds records itself under a stage-scoped key, and
            # that only works if the caller's RunSummary actually reaches here.
            seen.append((sequence, optimization_cub_index, run_summary))
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
        run_summary = RunSummary()

        patched_cai, patched_tai, summaries = main_module.run_hotspot_avoidance(
            module_input=module_input,
            cds_nt_final_cai=cai,
            cds_nt_final_tai=tai,
            skipped_codons_num=0,
            run_summary=run_summary,
        )

        assert [sequence for sequence, _, _ in seen] == cai + tai, (
            "every candidate in both lists must be patched"
        )
        # The CAI list must be scored with codon_adaptation_index and the tAI list with
        # trna_adaptation_index - never max_codon_trna_adaptation_index. That value would
        # resolve to a "max_cai_tai_profile" attribute Organism does not have, so
        # build_dcub_codon_table's getattr would fall back to {} and hand the optimizer an
        # all-zero codon table with no error raised anywhere.
        assert [index for _, index, _ in seen] == (
            [models.ORFOptimizationCubIndex.codon_adaptation_index] * len(cai)
            + [models.ORFOptimizationCubIndex.trna_adaptation_index] * len(tai)
        )
        # The caller's RunSummary must reach the module itself, not a stand-in:
        # the codon-loss table records itself there under a stage-scoped key.
        assert [summary for _, _, summary in seen] == [run_summary] * (len(cai) + len(tai))
        assert patched_cai == [candidate.replace("AAA", "AAG", 1) for candidate in cai]
        assert patched_tai == [candidate.replace("AAA", "AAG", 1) for candidate in tai]
        # The summary lookup is keyed by the PATCHED sequence, so the winner can
        # be resolved after evaluation picks it.
        for patched in patched_cai + patched_tai:
            assert summaries[patched]["num_edits"] == 1

    def test_translation_preserving_failure_is_caught_by_validate_module_output(self, monkeypatch):
        from modules import main as main_module
        from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult

        def broken_run_module(sequence, module_input, optimization_cub_index, skipped_codons_num,
                              run_summary, compute_motifs=None):
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
                run_summary=RunSummary(),
            )


class TestMotifDetectionFlag:
    """The UI's Advanced Options switch reaches ESO's motif scan.

    `None` is the meaningful default, not `False`: it means "no opinion from
    the request", which leaves HOTSPOT_AVOIDANCE.COMPUTE_MOTIFS in
    configuration.yaml in charge for callers that never send the field.
    """

    def test_module_input_defaults_to_deferring_to_the_config(self):
        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG",
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
        )
        assert module_input.enable_motif_detection is None

    def test_user_input_defaults_to_deferring_to_the_config(self):
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
        assert user_input.enable_motif_detection is None

    def test_user_input_accepts_both_overrides(self):
        for requested in (True, False):
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
                enable_motif_detection=requested,
            )
            assert user_input.enable_motif_detection is requested

    def test_module_input_summary_reports_the_flag(self):
        module_input = models.ModuleInput(
            organisms=[],
            sequence="ATG",
            output_path="",
            tuning_parameter=0.5,
            clusters_count=1,
            enable_motif_detection=True,
        )
        assert module_input.summary["enable_motif_detection"] is True

    def test_user_input_module_forwards_the_flag(self):
        from modules.user_IO.user_input import UserInputModule

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
            enable_motif_detection=True,
        )
        module_input = UserInputModule(user_input=user_input, skipped_codons_num=0).run_module(
            run_summary=RunSummary(),
        )
        assert module_input.enable_motif_detection is True

    def test_flag_reaches_the_hotspot_module(self, monkeypatch):
        from modules import main as main_module
        from modules.hotspot_avoidance.hotspot_avoidance_main import HotspotPatchResult

        seen = []

        def fake_run_module(sequence, module_input, optimization_cub_index, skipped_codons_num,
                            run_summary, compute_motifs=None):
            seen.append(compute_motifs)
            return HotspotPatchResult(
                sequence_before=sequence,
                sequence_after=sequence,
                num_edits=0,
                detected_sites={"recombination": 0, "slippage": 0, "motifs": 0},
                warnings=[],
            )

        monkeypatch.setattr(
            main_module.HotspotAvoidanceModule, "run_module", staticmethod(fake_run_module)
        )

        for requested in (True, False, None):
            seen.clear()
            module_input = models.ModuleInput(
                organisms=[],
                sequence="ATG" * 10,
                output_path="",
                tuning_parameter=0.5,
                clusters_count=1,
                enable_hotspot_avoidance=True,
                enable_motif_detection=requested,
            )
            main_module.run_hotspot_avoidance(
                module_input=module_input,
                cds_nt_final_cai=["AAA" * 10],
                cds_nt_final_tai=["CCC" * 10],
                skipped_codons_num=0,
                run_summary=RunSummary(),
            )
            assert seen == [requested, requested], (
                "both the CAI and the tAI candidate lists must carry the request's choice"
            )
