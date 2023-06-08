from pipline_steps_negative.duplex_step_negative import duplex as duplex_negative
from pipeline_steps.duplex_step import duplex as duplex_positive
from pipline_steps_negative.rna_site_insertion_negative import get_site_from_extended_site
from pipline_steps_negative.normalization_final_step_negative import finalize
from pipeline_steps.feature_extraction import feature_extraction
from consts.global_consts import ROOT_PATH,NEGATIVE_DATA_PATH, GENERATE_DATA_PATH
from utils.utilsfile import read_csv, to_csv
from pathlib import Path


# step 1- formalization all the dataset for the same format

def full_pipline(name_of_method: str, name_of_file: str, sign:int):
    # step 2- extract duplex of the interaction by VieannaDuplex
    # only for mock mirna
    name_of_file_primary = name_of_method + "_" + name_of_file.split('_negative')[0] + "_negative"

    # name_of_file_primary = name_of_method + "_" + name_of_file
    fout_primary = NEGATIVE_DATA_PATH / name_of_method

    name_of_file = name_of_file + ".csv"
    fin = GENERATE_DATA_PATH / name_of_method / name_of_file

    name_of_file = name_of_file_primary + "_duplex.csv"
    fout= fout_primary / name_of_file

    # positive interactions
    # if sign == 1:
    #     print("###############Duplex POSITIVE#############")
    #     duplex_positive('ViennaDuplex', fin, fout)
    # else:
    #     # negative interactions
    #     print("###############Duplex NEGATIVE#############")
    #     duplex_negative('ViennaDuplex', fin, fout)

    # step 3- extract the site and his coordination's
    fin = fout
    print("###############Site#############")

    name_of_file = name_of_file_primary + "_site.csv"
    fout = fout_primary /  name_of_file
    get_site_from_extended_site(fin, fout)

    print("###############Normaliztion#############")

    # step 4- normalization of the dataframe
    fin = fout
    name_of_file = name_of_file_primary + "_normalization.csv"
    fout = fout_primary / name_of_file
    finalize(fin, fout)

    # step 5- extract features
    print("###############extract features#############")

    fin = fout
    name_of_file = name_of_file_primary + "_features.csv"
    fout = fout_primary / name_of_file
    name_of_file_ch = name_of_file_primary + "_check_before_filter.csv"
    fout_check = fout_primary / name_of_file_ch
    feature_extraction(fin, fout)


def check_pipline():
    pos_file_name = "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv"
    pos_df = read_csv(pos_file_name)
    pos_df.rename(columns={"sequence": "full_mrna", "site": 'site_chimiri'}, inplace=True)
    col_list = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'site_chimiri', 'region', 'valid_row',
                'full_mrna', 'Gene_ID']
    pos_df = pos_df[col_list]
    to_csv(pos_df, Path(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/darnell_human_ViennaDuplex_features.csv/darnell_human_ViennaDuplex_check.csv"))
    # full_pipline("validation", "darnell_human_ViennaDuplex_check")

    pos = "/sise/home/efrco/efrco-master/data/negative_interactions/validation/validation_darnell_human_ViennaDuplex_negative_features.csv"
    open = read_csv(pos)
    open.drop(columns=['Seed_match_noncanonical', 'Seed_match_canonical'], inplace=True)
    to_csv(open, pos)

    ########################model###############################
    # from Classifier.train_test_split_stratify import split_train_test
    # from Classifier.ClassifierWithGridSearch import main_primary
    # from Classifier.result_test import different_results_summary
    # split_train_test()
    # main_primary()
    # different_results_summary(method_split="validation", model_dir="models_validation")


check_pipline()