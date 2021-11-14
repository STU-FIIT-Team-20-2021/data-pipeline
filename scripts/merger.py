"""Manages the merging function of the data pipeline."""
import pandas as pd


def merge() -> None:
    """
    Merges multiple .csv files.
    :return: None
    """
    file_pubchem = "..data/crawler_output/pubchem.csv"
    file_swiss = "..data/crawler_output/swiss.csv"
    out = "data/merger_output/merged.csv"

    df_swiss = pd.read_csv(file_swiss)
    df_pubchem = pd.read_csv(file_pubchem, index_col=0)

    df_concat = pd.concat([df_pubchem, df_swiss], axis=1)

    df_concat.to_csv(out, index=False, header=True)


if __name__ == "__main__":
    merge()
