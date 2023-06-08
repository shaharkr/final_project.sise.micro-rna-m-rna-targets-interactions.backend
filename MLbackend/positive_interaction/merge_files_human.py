import sys
from typing import Tuple
import pandas as pd
from numpy import int64
from pandas import DataFrame
from pathlib import Path
from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, ROOT_PATH, MERGE_DATA, POSITIVE_PATH_DATA
from consts.mirna_utils import MIRBASE_FILE
from utils import to_csv
from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, get_wrapper, to_csv, read_csv

# Q CLASH PRE-PROCESSING

def fill_full_name_mirRNA(file_name):

    frames = []
    for i in range(0, 5):
        new_file_name = str(file_name) + str(i) + ".csv"
        path_file_name = POSITIVE_PATH_DATA / new_file_name
        df: DataFrame = read_csv(path_file_name)
        df.reset_index(inplace=True)
        print(new_file_name, "  " , df.shape)
        df = filter(df)
        df = mirna_sequences(df)


        path = MERGE_DATA / "positive_interactions_merge/temp" / new_file_name
        df.drop(labels=['level_0'], axis=1, inplace=True)

        to_csv(df, path=path)


def mirna_sequences(df: DataFrame):

    df_origin = df
    file_name = MIRBASE_FILE
    logger.info(f"Reading file {file_name}")
    df_mirna = pd.read_csv(file_name)
    df_mirna_hsa = df_mirna[df_mirna['miRNA ID'].str.contains('hsa')]

    # In cases where multiple mirna exist because the different release , we considered the
    # last release.
    df_mirna_hsa['version'] = pd.to_numeric(df_mirna_hsa['version'])

    df_mirna_hsa = df_mirna_hsa.loc[df_mirna_hsa.groupby(['miRNA ID'])["version"].idxmax()]


    # for each mirna we find the longest sequence of the mirna
    df = pd.merge(df, df_mirna_hsa, how='inner', on=['miRNA sequence'])
    df['version'] = pd.to_numeric(df['version'])
    print("before", df.shape[0])
    df = df.loc[df.groupby(['index'])["version"].idxmax()]
    print("after", df.shape[0])

    df = df.drop(labels=['version', 'prefix', 'Unnamed: 0', 'miRNA ID_x'], axis=1)
    df.insert(4, "miRNA ID", df['miRNA ID_y'])
    df.drop(labels=['miRNA ID_y'], axis=1, inplace=True)
    return df


def merge_q_clash(file_name) -> DataFrame:

    frames = []
    dir = MERGE_DATA / "positive_interactions_merge/temp"
    files = list(dir.glob('**/*.csv'))
    for dataset in files:
        df: DataFrame = read_csv(dataset)
        # print(df.shape[0])
        frames.append(df)
    result = pd.concat(frames)

    result = result[result['sequence'].apply(lambda x: len(x)) > 40]


    result.reset_index(drop=True, inplace=True)
    result.reset_index(inplace=True)
    result.drop(labels=['level_0'], axis=1, inplace=True)

    return result



################################ALL-FIELS-PROCESSING#################################################

def merge(file_name) -> DataFrame:

    frames = []
    for i in range(0, 5):
        new_file_name = str(file_name)+str(i)+".csv"
        print(new_file_name)
        path_file_name = POSITIVE_PATH_DATA / new_file_name
        df: DataFrame = read_csv(path_file_name)
        print(df.shape[0])
        frames.append(df)
    result = pd.concat(frames)

    return result


def save(df: DataFrame, file_name: str):
    full_path = MERGE_DATA / "positive_interactions_merge" / file_name
    to_csv(df, full_path)


def filter(result) -> DataFrame:

    # # my append - 25/01
    # result = result[result['region count'] == 1]

    df = result[result['region'].str.contains('3utr')]
    df = df.rename(columns={"miRNA ID": "miRNAID"})

    # remove null miRNAID
    df = df[~df.miRNAID.isnull()]
    df = df.rename(columns={'miRNAID': 'miRNA ID'})

    #df = df[df['miRNA ID'].str.contains('hsa')]

    df = df[df['valid_row'] == True]
    print('3utr: ', "  ", df.shape)

    df = df[(df["Seed_match_canonical"] == 'True') | (df["Seed_match_noncanonical"] == 'True')]
    print('canon: ', "  ", df.shape)

    df = df[df['sequence'].apply(lambda x: len(x)) > 40]

    # df = df[df['miRNA sequence'].apply(lambda x: len(x)) > 21]
    df.drop(columns=['miRNAMatchPosition_21','miRNAMatchPosition_22'], inplace=True)

    df = df[~df.isna().any(axis=1)]

    df.reset_index(drop=True, inplace=True)
    df.reset_index(inplace=True)
    return df


def processing():

    # Dataset one
    # result_h1 = merge("human_mapping_ViennaDuplex_features")
    # result_h1 = filter(result_h1)
    # print("H1 after filter:", result_h1.shape)
    # save(result_h1, "human_mapping_ViennaDuplex_features.csv")
    #
    #
    # # Dataset two
    # result_h2 = merge("unambiguous_human_ViennaDuplex_features")
    # result_h2 = filter(result_h2)
    # print("H2 after filter:", result_h2.shape)
    # save(result_h2, "unambiguous_human_ViennaDuplex_features.csv")
    #
    # # # Data set three
    result_h3 = merge("darnell_human_ViennaDuplex_features")
    result_h3 = filter(result_h3)
    print("H3 after filter:", result_h3.shape)
    save(result_h3, "darnell_human_ViennaDuplex_features.csv")

    # Dataset four
    # fill_full_name_mirRNA("qclash_melanoma_human_ViennaDuplex_features")
    # result_h4 = merge_q_clash("qclash_melanoma_human_ViennaDuplex_features")
    # save(result_h4, "qclash_melanoma_human_ViennaDuplex_features.csv")

processing()



