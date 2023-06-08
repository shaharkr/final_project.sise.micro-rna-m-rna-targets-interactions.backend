from pathlib import Path
from datetime import datetime
from utils.logger import logger
from consts.mirna_utils import MIRBASE_FILE
from consts.global_consts import DUPLEX_DICT
from utils.utilsfile import read_csv, to_csv
from consts.global_consts import ROOT_PATH, DATA_PATH,CLIP_PATH_DATA,NEGATIVE_DATA_PATH, GENERATE_DATA_PATH
import pandas as pd
import numpy as np
from utils.utilsfile import *
from consts.global_consts import ROOT_PATH, DATA_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from duplex.ViennaDuplex import ViennaDuplex
from duplex.ViennaDuplex import *
import random
from features.SeedFeatures import *
# import MirBaseUtils.mirBaseUtils as MBU
import mirna_utils.mirbase as MBU
from multiprocessing import Process
from consts.global_consts import CONFIG
from duplex.Duplex import Duplex
import os


dir_after_filter = "split_after_filter"

def get_end(sequence, end, extra_chars):
    return min(len(sequence), end + extra_chars)

def sites_list_clash():
    extra_chars = 10
    positive_interaction_dir = MERGE_DATA / "positive_interactions_new/featuers_step/"
    df_sites = []
    # Gene_ID ---> number Gen + number Transcript
    # ID ---> number Gen
    for file in positive_interaction_dir.glob("*csv*"):
        logger.info(f"Reading file {file}")
        df_positive_interaction = read_csv(file)

        df_positive_interaction['site_extend'] = df_positive_interaction.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                                  "sequence", "start", "end", strand="+",
                                                  extra_chars=extra_chars), axis=1)

        df_positive_interaction['start'] = df_positive_interaction['start'].apply(lambda x: max(0, int(x) - 1 - extra_chars))
        df_positive_interaction['end'] = df_positive_interaction.apply(
            func=get_wrapper(get_end, "sequence", "end", extra_chars=extra_chars), axis=1)
        df_positive_interaction.rename(columns={'site_extend': 'site_old'}, inplace=True)
        df_sites.append(df_positive_interaction)

    df_merge = pd.concat(df_sites)
    df_merge = df_merge[['Gene_ID',  'site_old', 'start', 'end']]
    return df_merge

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
    duplex = RNA.duplexfold(mir, mrna)
    MEF_duplex = duplex.energy
    site = dp.site[::-1]

    return canonic_seed, non_canonic_seed, dp.interaction_count, MEF_duplex, site


def generate_negative_seq(full_mrna, orig_mirna, num_of_tries=10000):
    canonic_seed, non_canonic_seed, num_of_pairs, MEF_duplex, site = valid_negative_seq(orig_mirna, full_mrna)
    cond1 = canonic_seed
    cond2 = non_canonic_seed
    if cond1 or cond2:
        properties = {
            "mock_mirna": orig_mirna,
            "full_mrna": full_mrna,
            "canonic_seed": canonic_seed,
            "non_canonic_seed": non_canonic_seed,
            "num_of_pairs": num_of_pairs,
            "MEF_duplex": MEF_duplex,
            "site": site
        }
        return True, properties
    return False, {}


def generate_interactions(site, Gene_ID, full_mrna, df_mirna):
    valid = False
    num_tries = 0

    while not valid:
        random_numnber = random.randint(0, 1000)
        random_mirna = df_mirna.sample(n=1, random_state=random_numnber)
        num_tries += 1
        if num_tries > 100000:
            break
        random_mirna.reset_index(drop=True, inplace=True)
        random_mirna.reset_index(inplace=True)
        for index, row in random_mirna.iterrows():
            valid, properties = generate_negative_seq(site, row['sequence'])

            new_row = pd.Series()
            new_row['paper name'] = 'mirTarget'
            new_row['organism'] = 'Human'
            new_row['miRNA ID'] = row['miRNA ID']
            new_row['miRNA sequence'] = row['sequence']
            new_row['Gene_ID'] = Gene_ID
            new_row['full_mrna'] = full_mrna

    if valid:
        return new_row, properties['MEF_duplex'], valid
    else:
        return {}, {}, valid



def sub_insert_NNN(full_mrna, start, end, site):
    start = int(float(start))
    end = int(float(end))
    while start != end + 1:
        full_mrna = full_mrna[:start] + "N" + full_mrna[start + 1:]
        start += 1
    return full_mrna


def complete_site_chars(start, end):
    len_site = end - start
    number_chars_add_one_side = 0
    if len_site < 75:
        number_chars_add = 75 - len_site
        number_chars_add_one_side = round(number_chars_add / 2)

    return number_chars_add_one_side


def worker(df_mrna, df_sites, df_mirna):

    # group the interaction by the gene (find all the site of the same gene in the clip data)
    gruop_df = df_mrna.groupby(['Gene_ID'])
    neg_df = pd.DataFrame()
    count = 0
    for group_name, sub_group in gruop_df:
        mrna_cut = sub_group.iloc[0]["full_mrna"]
        Gene_ID = sub_group.iloc[0]["Gene_ID"]
        full_mrna = sub_group.iloc[0]["full_mrna"]
        # all the sites from the clash + CLIP
        df_sites_gene_id = df_sites[df_sites['Gene_ID'] == Gene_ID]

        # masking clip interactions
        # mask the original site of the interaction by the clip data
        for row_index, row in sub_group.iterrows():
            mrna_cut = sub_insert_NNN(mrna_cut, row["start"], row["end"], row['site_old'])

        # masking clash interactions
        for row_index, row in df_sites_gene_id.iterrows():
            mrna_cut = sub_insert_NNN(mrna_cut, row["start"], row["end"], row['site_old'])

        cut_mrna = mrna_cut
        size_param = 40
        pervious_MEF_duplex = float('inf')
        best_row = pd.Series()
        for window in range(0, len(cut_mrna) + size_param, size_param):
            current_row = []
            sub_mrna = cut_mrna[window:window+75]
            cd= len(cut_mrna)
            start_site = window
            end_site = window+75
            if end_site >= len(cut_mrna):
                break
            if "N" in sub_mrna:
                count = count + 1
                continue

            else:
                number_char_complete = int(complete_site_chars(start_site, end_site))
                full_sub = get_subsequence_by_coordinates(full_mrna, start_site,end_site, "+", number_char_complete)

                current_row, new_MEF_duplex, valid = generate_interactions(full_sub, Gene_ID, full_mrna, df_mirna)
                if not valid:
                    continue
                current_row['start'] = start_site
                current_row['end'] = end_site
                current_row['site'] = full_sub
                if new_MEF_duplex <= pervious_MEF_duplex:
                    best_row = current_row
                    pervious_MEF_duplex = new_MEF_duplex

        if len(best_row) == 0:
            count += 1
            print("not found interaction for")
            continue

        neg_df = neg_df.append(best_row, ignore_index=True)

    return neg_df

def clear_dir_genes():

    path = "/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/split_by_gene/"
    # Get all csv files in the specified directory
    files = [f for f in os.listdir(path) if f.endswith('.csv')]

    # Iterate over the files
    for file in files:
        # Delete the file
        os.remove(os.path.join(path, file))

def split_file_by_gene():

    clear_dir_genes()
    clip_data_path = CLIP_PATH_DATA
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*mrna_clean.csv*"):
            if "clip_3" not in str(file):
                continue
            name_mrna_file = str(clip_dir.stem) + ".csv"
            mrna_df = read_csv(GENERATE_DATA_PATH / "clip_interaction" / name_mrna_file)
            gruop_df = mrna_df.groupby(['Gene_ID'])
            number_file = 0
            for group_name, sub_group in gruop_df:
                name_dir = "non_overlapping_sites_clip_data"
                name_file = str(number_file) + "_" + str(group_name) + ".csv"
                path_df = GENERATE_DATA_PATH / name_dir / "split_by_gene" / name_file
                to_csv(sub_group, path_df)
                number_file = number_file + 1


def run(start, to):
    clip_data_path = CLIP_PATH_DATA
    list_sites_clash = sites_list_clash()

    usecols = ['Gene_ID', 'site_old', 'start', 'end']
    # add the fragments from clip experiments
    mrna_df_full = pd.read_csv(GENERATE_DATA_PATH / "clip_interaction" / 'clip_3.csv', usecols=usecols)
    all_clip_clash_site = list_sites_clash.append(mrna_df_full)

    frames = []
    for clip_dir in clip_data_path.iterdir():
        if "clip_3" not in str(clip_dir.stem):
            continue
        # arrive to the correct dir
        for file in clip_dir.glob("*mrna_clean.csv*"):
            mirna_df = read_csv(clip_dir / "mirna.csv")
            name_dir = "non_overlapping_sites_clip_data"
            path_df = GENERATE_DATA_PATH / name_dir / "split_by_gene/"

            # pass on the start until to files

            for file in path_df.glob("*.csv"):
                # name_mrna_file = str(start) + ".csv"

                files = [f for f in os.listdir(path_df) if f.endswith('.csv')]

                # Find the file whose name starts with the input name
                matching_file = [f for f in files if f.startswith(str(start))][0]

                mrna_df = read_csv(path_df / matching_file)
                neg_df = worker(mrna_df, all_clip_clash_site, mirna_df)
                frames.append(neg_df)
                start = start + 1
                if start == to + 1:
                    break

        name_dir = "non_overlapping_sites_clip_data"
        name_file = str(to) + ".csv"
        path_df = GENERATE_DATA_PATH / name_dir / dir_after_filter / name_file
        result = pd.concat(frames)
        print(result.shape)
        result.reset_index(drop=True, inplace=True)
        result.reset_index(inplace=True)
        to_csv(result, path_df)


def drop_duplicate(result):
    # Gene_ID ---> number Gen
    # ID ---> number Gen + number Transcript

    print("size before filter grop mrna:", result.shape)
    df_g = result.groupby(['ID', 'site_new'])
    result = result.loc[df_g["sequence length"].idxmax()]
    print("size after filter grop mrna:", result.shape)

    return result

def combine_files():

    name_dir = "non_overlapping_sites_clip_data"
    data_dir =  GENERATE_DATA_PATH / name_dir / dir_after_filter
    frames = []
    count = 0
    files = list(data_dir.glob("*.csv*"))
    for file_gene in data_dir.glob("*.csv*"):
        # if "clash" not in str(file_gene.name):
        #     continue
        df = read_csv(file_gene)
        count = count + df.shape[0]

        frames.append(df)

    print("f")
    print("COUNT: ######################################", count)
    result = pd.concat(frames)
    print(result.shape)
    result.reset_index(drop=True, inplace=True)

    result.reset_index(inplace=True)
    result.drop(columns=['level_0'], inplace=True)

    # result = drop_duplicate(result)
    name = "darnell_human_ViennaDuplex_features_negative" + ".csv"
    path_df = GENERATE_DATA_PATH / name_dir / name
    # result = result.reset_index(drop=True)
    # result.drop(columns=['level_0'], inplace=True)
    result.reset_index(drop=True, inplace=True)

    print('f')
    to_csv(result, path_df)



def main_run_clip():
    run(start=0, to=100)

    run(start=101, to=300)

    run(start=301, to=600)

    run(start=601, to=900)

    run(start=901, to=1000)

    run(start=1001, to=1200)

    run(start=1201, to=1400)

    run(start=1401, to=1600)

    run(start=1601, to=1815)

# split_file_by_gene()
# main_run_clip()
# combine_files()
# maybe we need to filter the interactions- because for each gene we have a lot of candidate
# check that for each gen we have only one candidate
# d = read_csv("/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/darnell_human_ViennaDuplex_features_negative.csv")
# print(d.shape)

def find_files_clash(name, path):
    results = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if name in file:
                results.append(file.split("_")[1])
    return results


def find_files_clip(name, path):
    results = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if name in file:
                results.append(file.split("_")[0])
    return results


def find_gene_CLIP_not_found_interactions_clip():
        clip_data_path = CLIP_PATH_DATA
        for clip_dir in clip_data_path.iterdir():
            for file in clip_dir.glob("*mrna_clean.csv*"):
                if "clip_3" not in str(file):
                    continue
                name_mrna_file = str(clip_dir.stem) + ".csv"
                mrna_df = read_csv(GENERATE_DATA_PATH / "clip_interaction" / name_mrna_file)
                k_clip_all= set(list(mrna_df['Gene_ID']))


        interaction_find = read_csv("/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/non_overlapping_sites_clip_data_2.csv")
        k_clip_sub = set(list(interaction_find['Gene_ID']))

        #list of gene that need rerun
        x = set(list(k_clip_all.difference(k_clip_sub)))
        all_files = []
        for gene in x:
            name = gene
            path = "/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/split_by_gene/"
            all_files.extend(find_files_clip(name, path))
        print(len(x))
        print(len(all_files))
        return all_files
#
# list_rerun = find_gene_CLIP_not_found_interactions_clip()
# print(list_rerun)


def clsah_split_gene():

    # clean directory
    path = "/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/split_by_gene_clash/"
    # Get all csv files in the specified directory
    files = [f for f in os.listdir(path) if f.endswith('.csv')]

    # Iterate over the files
    for file in files:
        # Delete the file
        os.remove(os.path.join(path, file))

    # generate file
    clip_data_path = CLIP_PATH_DATA
    clash_darnell = read_csv(
        "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")
    c = list(clash_darnell['Gene_ID'])
    # clash gene
    k_clash = set(c)
    for clip_dir in clip_data_path.iterdir():
        for file in clip_dir.glob("*mrna_clean.csv*"):
            if "clip_3" not in str(file):
                continue
            name_mrna_file = str(clip_dir.stem) + ".csv"
            mrna_df = read_csv(GENERATE_DATA_PATH / "clip_interaction" / name_mrna_file)
            # clip gen
            k_clip= set(list(mrna_df['Gene_ID']))

            clash_darnell_group = clash_darnell.groupby(['Gene_ID'])
            gene_clash = k_clash.difference(k_clip)
            gene_clip = k_clip.difference(k_clash)
            print("##########################:   ", len(gene_clash))
            number_file = 0
            for group_name, sub_group in clash_darnell_group:
                if group_name in gene_clash:
                    name_dir = "non_overlapping_sites_clip_data"
                    name_file = "clash_" + str(number_file) + "_" + str(group_name) + ".csv"
                    path_df = GENERATE_DATA_PATH / name_dir / "split_by_gene_clash" / name_file
                    to_csv(sub_group, path_df)
                    number_file = number_file + 1
            print(number_file)


def run_clash(start, to):
    clip_data_path = CLIP_PATH_DATA
    list_sites_clash = sites_list_clash()

    usecols = ['Gene_ID', 'site_old', 'start', 'end']
    mrna_df_full = pd.read_csv(GENERATE_DATA_PATH / "clip_interaction" / 'clip_3.csv', usecols=usecols)
    all_clip_clash_site = list_sites_clash.append(mrna_df_full)

    frames = []
    for clip_dir in clip_data_path.iterdir():
        if "clip_3" not in str(clip_dir.stem):
            continue
        # arrive to the correct dir
        for file in clip_dir.glob("*mrna_clean.csv*"):
            mirna_df = read_csv(clip_dir / "mirna.csv")
            name_dir = "non_overlapping_sites_clip_data"
            path_df = GENERATE_DATA_PATH / name_dir / "split_by_gene_clash/"

            # pass on the start until to files

            for file in path_df.glob("*.csv"):

                files = [f for f in os.listdir(path_df) if f.endswith('.csv')]

                # Find the file whose name starts with the input name
                start_name = "clash_" + str(start)
                matching_file = [f for f in files if f.startswith(start_name)][0]

                mrna_df = read_csv(path_df / matching_file)
                mrna_df.rename(columns={"sequence": "full_mrna", 'site': 'site_old'}, inplace=True )
                neg_df = worker(mrna_df, all_clip_clash_site, mirna_df)
                frames.append(neg_df)
                start = start + 1
                if start == to + 1:
                    break

        name_dir = "non_overlapping_sites_clip_data"
        name_file = "clash_duplicate" + str(to) + ".csv"
        path_df = GENERATE_DATA_PATH / name_dir / dir_after_filter / name_file
        result = pd.concat(frames)
        print(result.shape)
        result.reset_index(drop=True, inplace=True)
        result.reset_index(inplace=True)
        to_csv(result, path_df)


def main_run_clash():

    # run_clash(start=0, to=500)
    # run_clash(start=501, to=1000)
    # run_clash(start=1001, to=1500)
    run_clash(start=1501, to=1688)



# clsah_split_gene()
# main_run_clash()

# to append more interactions
def find_gene_apper_more_than_one_time():
        clash_darnell = read_csv(
            "/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_new/featuers_step/darnell_human_ViennaDuplex_features.csv")

        result = clash_darnell['Gene_ID'].value_counts()
        result = result[result > 2].index.tolist()
        all_files_twice = []
        for gene in result:
            name = gene
            path = "/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/split_by_gene_clash/"
            all_files_twice.extend(find_files_clash(name, path))

        find_interaction = read_csv("/sise/home/efrco/efrco-master/generate_interactions/non_overlapping_sites_clip_data/darnell_human_ViennaDuplex_features_negative.csv")
        number_to_choose = clash_darnell.shape[0] - find_interaction.shape[0]
        # return all_files_twice
        return all_files_twice, number_to_choose


def main_run_twice():
    all_files_twice, number_to_choose = find_gene_apper_more_than_one_time()
    for i in range(0, number_to_choose):
        chosen_number = random.choice(all_files_twice)
        all_files_twice.remove(chosen_number)
        run_clash(start=int(chosen_number), to=int(chosen_number))

# after that we need to check if the duplicate interaction are same
# main_run_twice()