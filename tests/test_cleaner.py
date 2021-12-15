"""Unit tests for cleaner.py"""
import pytest

from scripts.cleaner import *


class LoadFiles:
    instances = ["../data/merger_output/merged.csv",
                 "../data/crawler_output/swiss.csv",
                 "../data/crawler_output/pubchem.csv"]


@pytest.fixture(params=LoadFiles.instances)
def file_load(request):
    return request.param


def test_load_data(file_load):
    """
    Unit test. Checks load_data function.
    Scenario: checking if loading the file is ok.
    Assert: check if loaded data exists, if the file path is not string or file does not exist
    assert fails
    """
    if not isinstance(file_load, str):
        assert False
    try:
        df = load_data(file_load)
        assert df is not None
    except FileNotFoundError:
        assert False


def test_set_proper_data_type():
    """
    Unit test. Checks setting data type.
    Scenario: identifying whether the dataframe has changed its data types.
    Assert: dtypes comparison dataframe should not be empty.
    """
    df = load_data("../data/merger_output/merged.csv")
    df2 = set_proper_data_type(df)
    df_compare = df.dtypes.compare(df2.dtypes)
    assert not df_compare.empty
