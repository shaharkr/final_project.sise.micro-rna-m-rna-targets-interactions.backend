import pandas as pd
from pandas import DataFrame
from consts.global_consts import MERGE_DATA, ROOT_PATH, BIOMART_PATH
from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, to_csv
from consts.global_consts import ROOT_PATH, DATA_PATH,CLIP_PATH_DATA, GENERATE_DATA_PATH
from utils.utilsfile import read_csv, to_csv, get_wrapper
from consts.mirna_utils import MIRBASE_FILE
from consts.global_consts import DUPLEX_DICT
from duplex.ViennaDuplex import *
import random
from features.SeedFeatures import *
#import MirBaseUtils.mirBaseUtils as MBU
import mirna_utils.mirbase as MBU
from multiprocessing import Process
from consts.global_consts import CONFIG
from duplex.Duplex import Duplex

def get_all_seed_family():
    file_name = MIRBASE_FILE
    df_mirna = pd.read_csv(file_name)
    df_mirna_hsa = df_mirna[df_mirna['miRNA ID'].str.contains('hsa')]
    df_mirna_hsa['version'] = pd.to_numeric(df_mirna_hsa['version'])
    df_mirna_hsa = df_mirna_hsa.loc[df_mirna_hsa.groupby(['miRNA ID'])["version"].idxmax()]
    df_mirna_hsa['miRNA_ID_perfix'] = df_mirna_hsa['miRNA ID'].apply(lambda name: name.split("-")[2])
    df_mirna_hsa['seed_family'] = df_mirna_hsa['miRNA sequence'].apply(lambda m: extract_seed_family(m, []))
    df_mirna_hsa = df_mirna_hsa['seed_family'].unique()
    return list(df_mirna_hsa)

def change_columns_names(df: DataFrame) -> DataFrame:
    df = df.rename(columns={"mirna": "miRNA ID", "Gene_ID": "x"})
    df = df.rename(columns={"ID": "Gene_ID"})
    return df.rename(columns={"x": "ID"})


def add_meta_data(df: DataFrame) -> DataFrame:
    paper_name = "chimiRic"
    df.insert(0, "paper name", paper_name)
    df.insert(0, "organism", "Human")
    df["key"] = df.reset_index().index
    return df

def extract_seed_family(m: str, list_seed_family) -> str:
    try:
        seed = m[1:7]
        if seed not in list_seed_family:
            print("PROBLEM")
        return seed
    except TypeError:
        return ""

# In this step we remove all the interactions that exists in clash interaction
def filter_clash_interaction(df: DataFrame):

    positive_interaction_dir = MERGE_DATA / "positive_interactions_new/featuers_step/"
    print("################number of rows before filter clash interactions:#####", df.shape[0])

    # Gene_ID ---> number Gen + number Transcript
    # ID ---> number Gen

    for file in positive_interaction_dir.glob("*csv*"):
        usecols = ['miRNA ID', 'Gene_ID', 'seed_family']
        logger.info(f"Reading file {file}")
        df_positive_interaction = pd.read_csv(file, usecols=usecols)

        # transform the format of geneName
        df_positive_interaction['ID'] = df_positive_interaction['Gene_ID'].apply(lambda x: x.split("|")[0])

        # Intersection to find negative interaction wiche exists in clash poitive interacitons
        intersected_df = pd.merge(df, df_positive_interaction, how='inner', on=['seed_family', 'ID'])
        uniqueValues = len(intersected_df['key'].unique())
        # print("number of actual interaction to remove: ", uniqueValues)
        # remove interactions that exists in both of the dataset
        # print("number of rows to remove:" + str(file.stem) + " " + str(intersected_df.shape[0]))
        # print("number of rows before:" + str(df.shape[0]))
        new = df[(~df.key.isin(intersected_df.key))]
        # print("number of rows after:" + str(new.shape[0]))
        df = new

    print("################number of rows after filter clash interactions:##########", df.shape[0])

    return df


def merge_dataframes(mrna_file, mirna_file):
    mrna_file['key'] = 1
    mirna_file['key'] = 1
    # mrna_file = mrna_file.iloc[:3]
    # mirna_file = mirna_file.iloc[:3]
    df = mrna_file.merge(mirna_file, on='key').drop('key', axis=1)
    return df

def clean_dataframe(df):
    df = df.rename(columns={"site": "site_old"})
    df = df.drop(columns=['index_x', 'index_y', 'miRNA original', 'sequence_original', 'different', 'read_ID', 'key'], axis=1)
    df = df.rename(columns={"site_new": "site", "geneID": "Gene_ID" , 'sequence': 'miRNA sequence'})
    df["key"] = df.reset_index().index

    return df

def valid_negative_seq(mir, mrna):

    duplex_cls: Duplex = DUPLEX_DICT['ViennaDuplex']
    logger.info(f"{ViennaDuplex} do_duplex")
    dp = duplex_cls.fromChimera(mir, mrna)
    try:
        canonic_seed = dp.canonical_seed
        non_canonic_seed = dp.noncanonical_seed

    except SeedException:
        canonic_seed = False
        non_canonic_seed = False

    # warning: number of pair was before to interactions count
    return canonic_seed, non_canonic_seed, dp.interaction_count

def valid_interactions(mirna_sequence, mrna_site_sequence):

    # this step check that the mock is vaild. The meaning is that the interaction is canonical
    # or noncanonical

    canonic_seed, non_canonic_seed, num_of_pairs = valid_negative_seq(mirna_sequence, mrna_site_sequence)
    cond1 = canonic_seed
    cond2 = non_canonic_seed
    if cond1 or cond2:
        return True
    return False

def generate_interaction(df_mrna,mirna_df,list_seed_family):

    merge_interaction = merge_dataframes(df_mrna, mirna_df)
    merge_interaction = change_columns_names(merge_interaction)
    merge_interaction = add_meta_data(merge_interaction)
    merge_interaction['seed_family'] = merge_interaction['sequence'].apply(lambda m:extract_seed_family(m, list_seed_family))

    # step-1 ---> filter clash interaction
    merge_interaction_non_clash = filter_clash_interaction(merge_interaction)

    # step-2 ---> filter interactions by conditions - canon or non canon
    # site_new--> the original site after extracted from the full mrna by coordinates
    merge_interaction_non_clash["valid"] = merge_interaction_non_clash.apply(func=get_wrapper(valid_interactions,
                                                       "sequence", "site_new"), axis=1)

    df_final = merge_interaction_non_clash[merge_interaction_non_clash['valid'] == True]
    return df_final


def genetate_interactions_from_files():
    # split the dataframe to n frame
    clip_data_path = CLIP_PATH_DATA
    list_seed_family = get_all_seed_family()
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*mrna.csv*"):
            if "clip_3" not in str(file):
                continue
            mirna_df = read_csv(clip_dir / "mirna.csv")
            mrna_df = read_csv(clip_dir / "mrna_clean.csv")
            interactions_df = generate_interaction(mrna_df, mirna_df, list_seed_family)
            interactions_df = clean_dataframe(interactions_df)
            name_dir = "clip_interaction"
            name_file = str(clip_dir.stem) + ".csv"

            path_df = GENERATE_DATA_PATH / name_dir / name_file
            to_csv(interactions_df, path_df)

genetate_interactions_from_files()

