from consts.global_consts import ROOT_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
import pandas as pd
# from utils.utilsfile import read_csv, to_csv
from utilsfile import read_csv, to_csv

def tarBase_experminents():
    dir = NEGATIVE_DATA_PATH / "tarBase"
    name_file = "tarBase_human_negative_features.csv"
    df_tarBase = read_csv(dir / name_file)

    df_filter_Liver= df_tarBase[df_tarBase['tissue'] == "Liver"]
    df_filter_microArray= df_tarBase[df_tarBase['method'] == "Microarrays"]

    name_Liver = "tarBase_Liver_human_negative_features.csv"
    name_microArray = "tarBase_microArray_human_negative_features.csv"

    to_csv(df_filter_Liver, dir/name_Liver)
    to_csv(df_filter_microArray, dir/name_microArray)

# tarBase_experminents()


def clip_experminents():
    dir = NEGATIVE_DATA_PATH / "clip_interaction"
    name_file = "clip_interaction_clip_3_negative_features.csv"
    df_clip = read_csv(dir / name_file)
    df_filter_tail = df_clip[df_clip['not_match_site'] < 5]
    name_filter_tail = "clip_filter_tail_interaction_clip_3_negative_features.csv"
    to_csv(df_filter_tail, dir/name_filter_tail)


# clip_experminents()
