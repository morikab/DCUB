"""Regression test for app/modules/user_IO/_codonbias_compat.py.

Locks in the pandas>=3.0 fix for codon-bias's process_GtRNAdb_table: builds
a synthetic table with the same column dtypes pandas>=3.0's read_html
actually produces (pandas "string" dtype for text columns, int64 for a
numeric column mixed in, matching the real structure of a scraped GtRNAdb
page - confirmed by scraping a real page directly), then asserts the
ORIGINAL codon-bias implementation still reproduces the bug it was written
against (so this test would fail loudly if codon-bias ever fixes this
upstream, signaling the shim can be deleted), and that the patched version
- already active via importing modules.user_IO._codonbias_compat, exactly
as input_functions.py does - correctly extracts anti-codon/GCN pairs from
it.
"""
import codonbias.utils as cb_utils
import pandas as pd
import pytest

from modules.user_IO import _codonbias_compat

SYNTHETIC_GTRNADB_TABLE = pd.DataFrame({
    "Isotype": pd.array(["Ala", "Cys"], dtype="string"),
    "tRNA Count by Anticodon": pd.array(["AGC", "GGC 1"], dtype="string"),
    "tRNA Count by Anticodon.1": pd.array(["CGC", "TGC 5"], dtype="string"),
    "Total": pd.array([1, 2], dtype="int64"),
})


def test_original_codonbias_implementation_still_reproduces_the_bug():
    with pytest.raises(AttributeError, match="Can only use .str accessor"):
        _codonbias_compat._original_process_GtRNAdb_table(SYNTHETIC_GTRNADB_TABLE)


def test_patched_implementation_extracts_anti_codon_gcn_pairs():
    result = cb_utils.process_GtRNAdb_table(SYNTHETIC_GTRNADB_TABLE)
    pairs = dict(zip(result["anti_codon"], result["GCN"]))
    assert pairs == {"GGC": 1, "TGC": 5}
