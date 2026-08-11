"""Compatibility shim for codon-bias==0.3.1's process_GtRNAdb_table, which
breaks under pandas>=3.0.

Root cause: codon-bias filters an HTML table's columns via `table.dtypes ==
object` to find the text columns it needs to parse (anti-codon/count pairs
scraped from GtRNAdb). pandas>=3.0 changed `read_html`'s default dtype for
parsed string columns from the classic numpy `object` dtype to the newer
pandas `StringDtype` ("string") extension dtype, so that filter now silently
matches zero columns instead of raising. The resulting empty selection
cascades into a spuriously float64-typed empty "pair" column a few lines
later, which then crashes with `AttributeError: Can only use .str accessor
with string values, not floating` the moment `.str` is used on it. Confirmed
directly against live GtRNAdb data for both Escherichia coli and Bacillus
subtilis (the two organisms this codebase supports tAI weights for - see
calculate_tai_weights in input_functions.py). Also confirmed against
codon-bias==0.5.0 (the latest release as of this fix) - it has the exact
same `dtypes == object` logic, so upgrading does not fix this on its own.

This patches in a corrected column selector - `select_dtypes(include=
["object", "string"])`, which matches both the old and new pandas
string-dtype names - so tAI weight calculation keeps working.
`_original_process_GtRNAdb_table` is kept importable so a test can assert
the original still reproduces the bug this shim was written against;
delete this whole module (and its import in input_functions.py) once
codon-bias ships a real fix upstream.
"""
import codonbias.utils as _cb_utils
import pandas as pd

_original_process_GtRNAdb_table = _cb_utils.process_GtRNAdb_table


def _process_gtrnadb_table_pandas3_safe(table: pd.DataFrame) -> pd.DataFrame:
    df = table.select_dtypes(include=["object", "string"]).apply(
        lambda col: col.str.split(" ").str[-2:]
    )
    df = pd.DataFrame({"pair": df.values[df.apply(lambda col: col.str.len()).values == 2]})
    df["anti_codon"] = df["pair"].str[0]
    df["GCN"] = df["pair"].str[1].str.split("/").apply(lambda x: sum(map(int, x)))
    return df.drop(columns="pair")


_cb_utils.process_GtRNAdb_table = _process_gtrnadb_table_pandas3_safe
