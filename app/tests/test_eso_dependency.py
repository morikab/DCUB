"""ESO is a hard runtime dependency of the hotspot-avoidance module. It is an
unpublished sibling repo wired in as a path dependency, which is exactly the
kind of thing that silently breaks on a fresh clone - so assert it imports."""


def test_eso_public_api_is_importable():
    from eso import suspect_site_extractor
    from eso.optimize import optimization_engine

    assert callable(suspect_site_extractor)
    assert callable(optimization_engine)


def test_optimization_engine_accepts_the_parameters_this_integration_relies_on():
    import inspect

    from eso.optimize import optimization_engine

    parameters = inspect.signature(optimization_engine).parameters
    for name in (
        "custom_score_fn",
        "custom_score_minimize",
        "exclusion_regions",
        "orf_regions",
        "df_recombination",
        "df_slippage",
        "df_motifs",
        "mini_gc",
        "maxi_gc",
    ):
        assert name in parameters, f"eso.optimize.optimization_engine lost the {name!r} parameter"
