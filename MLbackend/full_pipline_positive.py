from utils.utilsfile import read_csv, to_csv
from pathlib import Path
from pipeline_steps.duplex_step import duplex as duplex_positive
from pipline_steps_negative.rna_site_insertion_negative import get_site_from_extended_site
from pipline_steps_negative.normalization_final_step_negative import finalize
from pipeline_steps.feature_extraction import feature_extraction
from consts.global_consts import MERGE_DATA,NEGATIVE_DATA_PATH, GENERATE_DATA_PATH, ROOT_PATH


# step 1- formalization all the dataset for the same format

def full_pipline(input_data):

    # positive interactions
    print("###############Duplex POSITIVE#############")
    duplex_positive_res = duplex_positive('ViennaDuplex', input_data)

    # step 3- extract the site and his coordination's
    print("###############Site#############")

    get_site_from_extended_site_res = get_site_from_extended_site(duplex_positive_res)

    print("###############Normaliztion#############")

    # step 4- normalization of the dataframe
    finalize_res = finalize(get_site_from_extended_site_res)

    # step 5- extract features
    print("###############extract features#############")

    feature_extraction_res = feature_extraction(finalize_res)
    return feature_extraction_res


# def generate_positive_interaction():
#
#     pos_dir_name = MERGE_DATA / "positive_interactions_merge"
#     for dataset_file in pos_dir_name.glob("*_features*"):
#         name_darnell = "darnell_human_ViennaDuplex_features"
#         name_data = str(dataset_file.stem)
#         if name_darnell != name_data:
#             continue
#         print(dataset_file)
#         pos_df = read_csv(dataset_file)
#         pos_df.rename(columns={"sequence": "full_mrna"}, inplace=True)
#
#         col_list = ['key', 'paper name', 'organism', 'miRNA ID', 'miRNA sequence', 'site', 'region','valid_row' , 'full_mrna', 'Gene_ID', 'region count']
#         pos_df = pos_df[col_list]
#         path = MERGE_DATA / "positive_interactions_new/data_without_featuers"
#         dataset_name = str(dataset_file.stem).split("_features.csv")[0].split("_features")[0]
#
#         name_file = path / (dataset_name + ".csv")
#         to_csv(pos_df, name_file)
#         # if str(dataset_name) == 'unambiguous_human_ViennaDuplex':
#         #     continue
#
#         print("full pipline for : ", dataset_file)
#         full_pipline(dataset_name)
#
#         # pos = MERGE_DATA / "positive_interactions_new/featuers_step" / (str(dataset_file.stem)+'.csv')
#         # open = read_csv(pos)
#         #
#         # # open.drop(columns=['Seed_match_noncanonical', 'Seed_match_canonical'], inplace=True)
#         # open = open[~open.isna().any(axis=1)]
#         # to_csv(open, pos)
#
# generate_positive_interaction()



