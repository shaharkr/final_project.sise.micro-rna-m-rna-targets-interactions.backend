from pathlib import Path
from typing import List
from functools import reduce
from consts.global_consts import NEGATIVE_DATA_PATH, ROOT_PATH, BIOMART_PATH
import pandas as pd
from pandas import DataFrame, Series
from consts.global_consts import SITE_EXTRA_CHARS

# from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, get_subsequence_by_coordinates_no_exception, \
    get_substring_index, get_wrapper, read_csv, to_csv

def extract_seed_family(m: str) -> str:
    try:
        return m[1:7]
    except TypeError:
        return ""



def finalize(df: DataFrame):

    def eta(x):
        try:
            return int(x) - SITE_EXTRA_CHARS
        except Exception:
            print(x)
            raise Exception()

    df['valid_row'] = True

    # logger.info("replace T with U")
    seq_cols = ['miRNA sequence', 'site', 'full_mrna', 'mrna_bulge', 'mir_bulge', 'mrna_inter', 'mir_inter' ]
    df[seq_cols] = df[seq_cols].replace(to_replace='T', value='U', regex=True)
    # logger.info("Add seed family")
    df["seed_family"] = df['miRNA sequence'].apply(extract_seed_family)

    # logger.info("Add valid/invalid flag")
    invalid_conditions = [pd.isna(df["miRNA sequence"]),
                          pd.isna(df["site"]),
                          df["miRNA sequence"].str.contains('X'),
                          df["miRNA sequence"].str.contains('N'),
                          df["site"].str.contains("N"),
                          df["site"].str.contains("Error"),
                          df["full_mrna"].str.contains('N'),
                          df["full_mrna"].str.contains('X'),
                          df["full_mrna"].str.contains("Error"),
                          df["full_mrna"].str.contains("None")]
    df["valid_row"] = ~reduce((lambda x, y: x | y), invalid_conditions)
    df["region"] = '3utr'
    df.reset_index(inplace=True)
    columns_name= list(df.columns)
    if 'level_0' in columns_name:
       df.drop(columns=['level_0'], inplace=True)

    return df




if __name__ == '__main__':
    pass
    # cli()
    # fin = NEGATIVE_DATA_PATH / "tarBase/duplex_rna_site.csv"
    # fout = NEGATIVE_DATA_PATH / "tarBase/tarBase_human_negative_normalization.csv"
    # #finalize(fin, fout)
    # df = pd.read_csv(Path(fout))
    # print(df)

