"""The UI posts hand-written literals straight into UserInput.

Nothing else checks the two sides agree: the frontend has no view of the
Python enums, and the backend never sees the TSX. A mismatch is invisible
until a user picks the affected option and FastAPI answers 422 - which is how
every Z-Score method shipped broken from the first Electron commit (273074d),
sending `zscore_bulk_diff` where the enum says `zscore_bulk_aa_diff`. Reading
the literals out of the source is ugly, but it is the only place the contract
can be checked at all.

Skipped when the UI sources are absent (a backend-only checkout, or a
packaged build) rather than failing there.
"""

import json
import re
from pathlib import Path

import pytest

from modules import models

UI_DIRECTORY = Path(__file__).resolve().parents[2] / "ui" / "DCUB"
ADVANCED_OPTIONS_PANEL = UI_DIRECTORY / "components" / "advanced-options-panel.tsx"
PAGE = UI_DIRECTORY / "app" / "page.tsx"

pytestmark = pytest.mark.skipif(
    not (ADVANCED_OPTIONS_PANEL.exists() and PAGE.exists()),
    reason="UI sources are not present in this checkout",
)


def _optimization_method_values() -> list:
    """The `optimizationMethods` dropdown's value strings, in source order."""
    source = ADVANCED_OPTIONS_PANEL.read_text(encoding="utf-8")
    block = re.search(
        r"const optimizationMethods = \[(.*?)\n  \]", source, re.DOTALL
    )
    assert block, "optimizationMethods array not found - has the panel been restructured?"
    values = re.findall(r'\{\s*value:\s*"([^"]+)"', block.group(1))
    assert values, "no dropdown values parsed"
    return values


def _cub_index_values() -> list:
    """The CUB index radio group's value strings."""
    source = ADVANCED_OPTIONS_PANEL.read_text(encoding="utf-8")
    return re.findall(r'<RadioGroupItem value="([^"]+)" id="(?:cai|tai)"', source)


def _request_payload_literals() -> dict:
    """The constant string literals in the `user_input_dict` payload."""
    source = PAGE.read_text(encoding="utf-8")
    block = re.search(r"user_input_dict: \{(.*?)\n        \},", source, re.DOTALL)
    assert block, "user_input_dict payload not found - has page.tsx been restructured?"
    return dict(re.findall(r'^\s*(\w+):\s*"([^"]*)",', block.group(1), re.MULTILINE))


#: Backend methods the dropdown deliberately does not offer. single_wanted
#: _organism scores against one wanted organism's profile alone, which sits
#: outside the wanted-vs-unwanted tradeoff the rest of the UI is built around.
#: Listed here so the coverage test below stays a real check rather than a
#: rubber stamp - if a new method lands in the enum, that test fails until it
#: is either offered or added here on purpose.
UNOFFERED_METHODS = {"single_wanted_organism"}


def test_the_dropdown_offers_every_backend_method_bar_documented_gaps():
    """All three weakest-link methods were missing from the dropdown while
    being fully implemented and working end to end (verified on mCherry against
    E. coli/B. subtilis, with and without hotspot avoidance)."""
    valid = {method.value for method in models.ORFOptimizationMethod}
    missing = valid - set(_optimization_method_values()) - UNOFFERED_METHODS
    assert not missing, f"backend methods the UI does not offer: {sorted(missing)}"


def test_every_optimization_method_option_is_a_real_enum_value():
    valid = {method.value for method in models.ORFOptimizationMethod}
    offered = _optimization_method_values()
    rejected = [value for value in offered if value not in valid]
    assert not rejected, (
        f"the UI offers optimization methods the backend rejects with a 422: "
        f"{rejected}. Valid values: {sorted(valid)}"
    )


def test_the_default_optimization_method_is_offered_and_valid():
    """store.ts's initial value and the panel's Reset to Defaults both name a
    method; a default the backend rejects breaks a first-time user's very
    first run."""
    store = (UI_DIRECTORY / "lib" / "store.ts").read_text(encoding="utf-8")
    default = re.search(r'optimizationMethod: "([^"]+)"', store)
    assert default, "no default optimizationMethod found in store.ts"
    assert default.group(1) in {method.value for method in models.ORFOptimizationMethod}
    assert default.group(1) in _optimization_method_values()


def test_every_cub_index_option_is_a_real_enum_value():
    valid = {index.value for index in models.ORFOptimizationCubIndex}
    offered = _cub_index_values()
    assert offered, "no CUB index options parsed"
    assert set(offered) <= valid, f"unknown CUB index values: {set(offered) - valid}"


def test_constant_payload_literals_match_their_enums():
    payload = _request_payload_literals()
    assert payload.get("initiation_optimization_method") in {
        method.value for method in models.InitiationOptimizationMethod
    }
    assert payload.get("evaluation_score") in {
        score.value for score in models.EvaluationScore
    }


def test_payload_keys_are_all_fields_userinput_accepts():
    """Pydantic ignores unknown keys by default, so a typo'd or renamed field
    does not 422 - it silently falls back to the model default. Catch it here
    instead."""
    source = PAGE.read_text(encoding="utf-8")
    block = re.search(r"user_input_dict: \{(.*?)\n        \},", source, re.DOTALL)
    sent = set(re.findall(r"^\s*(\w+):", block.group(1), re.MULTILINE))

    # pydantic v1 (see pyproject) - __fields__, and .alias is the field name
    # itself when no alias was declared.
    accepted = set()
    for name, field in models.UserInput.__fields__.items():
        accepted.add(name)
        accepted.add(field.alias)

    assert sent <= accepted, f"UserInput has no such field(s): {sorted(sent - accepted)}"


def test_persisted_method_renames_cover_every_stale_value():
    """The migration in store.ts must map every value the old build could have
    persisted; one left out keeps 422-ing for anyone who used that option."""
    store = (UI_DIRECTORY / "lib" / "store.ts").read_text(encoding="utf-8")
    block = re.search(
        r"const RENAMED_OPTIMIZATION_METHODS: Record<string, string> = \{(.*?)\n\}",
        store,
        re.DOTALL,
    )
    assert block, "RENAMED_OPTIMIZATION_METHODS not found in store.ts"
    renames = dict(re.findall(r"(\w+):\s*\"([^\"]+)\"", block.group(1)))
    valid = {method.value for method in models.ORFOptimizationMethod}

    assert renames, "no renames parsed"
    for stale, replacement in renames.items():
        assert stale not in valid, f"{stale} is a valid value; it needs no rename"
        assert replacement in valid, f"{stale} is migrated to unknown value {replacement}"
        assert replacement in _optimization_method_values(), (
            f"{stale} is migrated to {replacement}, which the dropdown no longer offers"
        )


def test_the_full_ui_payload_validates_against_userinput():
    """End to end on the shape the UI actually posts, for every method it
    offers - the check the 422 came from."""
    payload_template = {
        "sequence_file_path": None,
        "sequence": "ATGGCTTGTGATGAACATATCAAGCTGAAT",
        "tuning_param": 0.5,
        "organisms": {},
        "clusters_count": 1,
        "orf_optimization_cub_index": "CAI",
        "initiation_optimization_method": "original",
        "output_path": "results/DCUB/1",
        "evaluation_score": "average_distance",
        "enable_hotspot_avoidance": True,
        "enable_motif_detection": False,
    }
    for method in _optimization_method_values():
        user_input = models.UserInput(
            **json.loads(json.dumps({**payload_template, "orf_optimization_method": method}))
        )
        assert user_input.orf_optimization_method.value == method
