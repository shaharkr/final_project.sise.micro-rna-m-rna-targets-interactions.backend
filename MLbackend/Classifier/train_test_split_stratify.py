import pandas as pd
import numpy as np
from pandas import Series
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import mutual_info_classif
from sklearn.utils import shuffle
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Process
from consts.global_consts import DATA_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
# from utils.utilsfile import read_csv, to_csv
from utilsfile import read_csv, to_csv

def stratify_train_test_split(df, test_size , random_state):

    # Change: all the unique miRNA were put in the test set
    df.rename(columns={"miRNA ID": "microRNA_name"}, inplace=True)


    print("###############################################################")
    print(df.groupby("microRNA_name").ngroups)

    uniques_mirna = df[~df.duplicated(subset=['microRNA_name'], keep=False)]
    non_uniques_mirna = df[df.duplicated(subset=['microRNA_name'], keep=False)]

    # data for print
    non_uniques_mirna_group_by = non_uniques_mirna.groupby("microRNA_name")

    # dealing with the non_uniques_mirna
    non_uniques_train, non_uniques_test = train_test_split(non_uniques_mirna, test_size=test_size,
                                                           random_state=random_state,
                                                           stratify=non_uniques_mirna["microRNA_name"])

    train = pd.concat([non_uniques_train])
    test = pd.concat([non_uniques_test, uniques_mirna])
    train.rename(columns={"microRNA_name": "miRNA ID"}, inplace=True)
    test.rename(columns={"microRNA_name": "miRNA ID"}, inplace=True)

    return train, test


def stratify_train_test_split_worker(pos_train, pos_test, input_file_negative, output_dir, test_size, random_state):

    print(f"working on file: neg={input_file_negative}")
    neg = read_csv(input_file_negative)
    neg.insert(0, "Label", 0)
    print("neg", input_file_negative, neg.shape)

    df = set(pos_train.columns.values) - set(neg.columns.values)
    print(df, "CHANGE")

    # Both dataset must have the same columns
    col = [c for c in pos_train.columns if c in neg.columns]
    col_no = [c for c in pos_train.columns if c not in neg.columns]
    print("no col:", col_no)

    pos_train = pos_train[col]
    pos_test = pos_test[col]
    neg = neg[col]

    neg_train, neg_test = stratify_train_test_split(neg, test_size, random_state)

    # concat the pos & neg
    output = {
        "train": pos_train.append(neg_train, ignore_index=True),
        "test": pos_test.append(neg_test, ignore_index=True)
    }

    # save to csv
    dataset = str(input_file_negative.stem).split("_features_")[0]
    for key, d in output.items():
        out_file = f"{dataset}_{key}"
        out_file = out_file + "_stratify_method.csv"
        d = d.reindex(np.random.RandomState(seed=random_state).permutation(d.index))
        output_dir_train_test = output_dir / key / "stratify"
        to_csv(d, output_dir_train_test / out_file)



def test_size():

    print("###############TEST###################")
    dir = DATA_PATH_INTERACTIONS / "test/stratify"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*test*"):
        dataset = str(dataset_file.stem).split("_test")[0]
        for suffix in ["test_random_method", "test_stratify_method"]:
            if suffix == "test_random_method":
                dir = DATA_PATH_INTERACTIONS / "test/random"
            else:
                dir = DATA_PATH_INTERACTIONS / "test/stratify"
            d = read_csv(dir / f"{dataset}_{suffix}.csv")
            suffix = suffix.replace("test_", "")
            result.loc[dataset, suffix] = int(d.shape[0])
    result.sort_index(inplace=True)
    result = result.astype('int')
    print(result)
    print(60 * "*")
    # print(result.to_latex())


def train_size():

    print("###############Train###################")
    dir = DATA_PATH_INTERACTIONS / "train/stratify"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*train*"):
        dataset = str(dataset_file.stem).split("_train")[0]
        for suffix in ["train_random_method", "train_stratify_method"]:
            if suffix == "train_random_method":
                dir = DATA_PATH_INTERACTIONS / "train/random"
            else:
                dir = DATA_PATH_INTERACTIONS / "train/stratify"
            d = read_csv(dir / f"{dataset}_{suffix}.csv")
            suffix = suffix.replace("train_", "")
            result.loc[dataset, suffix] = int(d.shape[0])
    result.sort_index(inplace=True)
    result = result.astype('int')
    print(result)
    print(60 * "*")
    # print(result.to_latex())


# function that response to split with random
def split_train_test():
    train_test_dir = DATA_PATH_INTERACTIONS
    dir = NEGATIVE_DATA_PATH

    data_set_positive = MERGE_DATA / "positive_interactions_new/featuers_step" / "darnell_human_ViennaDuplex_features2.csv"

    # first we dived the positive interactions, the partition is const for all the train and test set.
    print(f"working on file: pos={data_set_positive}")
    pos = read_csv(data_set_positive)
    pos.insert(0, "Label", 1)
    print("pos", data_set_positive, pos.shape)
    pos_train, pos_test = stratify_train_test_split(pos, 0.2, 1*19)
    pos_train.rename(columns={"microRNA_name": "miRNA ID"}, inplace=True)
    pos_test.rename(columns={"microRNA_name": "miRNA ID"}, inplace=True)

    for method_dir in dir.iterdir():
        print(method_dir)
        for dataset_file_neg in method_dir.glob("*features*"):
           stratify_train_test_split_worker(pos_train, pos_test, dataset_file_neg, train_test_dir, 0.2, 1*19)


split_train_test()
# test_size()
# train_size()
# data_set_positive = MERGE_DATA / "darnell_human_ViennaDuplex_features.csv"
# pos = read_csv(data_set_positive)
# print(pos.shape[0])
# dir = NEGATIVE_DATA_PATH / "tarBase" / "tarBase_human_negative_features.csv"
# read_csv(dir)