"""Unit tests for merger.py"""
import pandas as pd

from scripts.merger import merge


def test_merge():
    """
    Unit test. Checks the merge function.
    Scenario: merging of pubchem and swiss csv files.
    Assert: equality of lists containing column names (merged.csv and result constructed here).
    """
    file_pubchem = "../data/crawler_output/pubchem.csv"
    file_swiss = "../data/crawler_output/swiss.csv"
    file_out = "../data/merger_output/merged.csv"
    cols_pubchem = pd.read_csv(file_pubchem, index_col=0, nrows=0).columns.tolist()
    cols_swiss = pd.read_csv(file_swiss, nrows=0).columns.tolist()
    merge()
    cols_out = pd.read_csv(file_out, nrows=0).columns.tolist()
    result = cols_pubchem
    result.extend(x+".1" if x in result else x for x in cols_swiss)
    assert result == cols_out
