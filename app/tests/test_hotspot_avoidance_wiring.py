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
