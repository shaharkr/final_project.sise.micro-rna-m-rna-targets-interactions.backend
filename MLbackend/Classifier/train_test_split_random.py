import pandas as pd
import numpy as np
from collections import Counter
from sklearn.datasets import make_classification
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import mutual_info_classif
from sklearn.utils import shuffle
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Process
from consts.global_consts import DATA_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from utilsfile import read_csv, to_csv

def random_train_test_split_worker(pos_train, pos_test, input_file_negative, output_dir, test_size, random_state):


    print(f"working on file: neg={input_file_negative}")
    neg = read_csv(input_file_negative)
    neg.insert(0, "Label", 0)

    # Both dataset must have the same columns
    col = [c for c in pos_train.columns if c in neg.columns]
    pos_train = pos_train[col]
    pos_test = pos_test[col]
    neg = neg[col]

    #pos_train, pos_test = train_test_split(pos, test_size=test_size, random_state=random_state)
    neg_train, neg_test = train_test_split(neg, test_size=test_size, random_state=random_state)

    # concat the pos & neg
    output = {
        "train": pos_train.append(neg_train, ignore_index=True),
        "test": pos_test.append(neg_test, ignore_index=True)
    }

    # save to csv
    dataset = str(input_file_negative.stem).split("_features_")[0]
    for key, d in output.items():
        out_file = f"{dataset}_{key}_random_method.csv"
        d = d.reindex(np.random.RandomState(seed=random_state).permutation(d.index))
        output_dir_train_test = output_dir / key / "random"
        to_csv(d, output_dir_train_test / out_file)


# function that response to split with random
def split_train_test():
    dir = NEGATIVE_DATA_PATH
    train_test_dir = DATA_PATH_INTERACTIONS
    data_set_positive = MERGE_DATA / "darnell_human_ViennaDuplex_features.csv"

    # first we dived the positive interactions, the partition is const for all the train and test set.
    print(f"working on file: pos={data_set_positive}")
    pos = read_csv(data_set_positive)
    pos.insert(0, "Label", 1)
    print("pos", data_set_positive, pos.shape)
    pos_train, pos_test = train_test_split(pos, test_size=0.2, random_state=19)

    for method_dir in dir.iterdir():
        print(method_dir)
        for dataset_file_neg in method_dir.glob("*features*"):
           random_train_test_split_worker(pos_train, pos_test, dataset_file_neg, train_test_dir, 0.2, 1*19)


def size():
    dir = NEGATIVE_DATA_PATH
    train_test_dir = DATA_PATH_INTERACTIONS
    result = pd.DataFrame()
    data_set_positive = MERGE_DATA / "darnell_human_ViennaDuplex_features.csv"

    # first we dived the positive interactions, the partition is const for all the train and test set.
    pos = read_csv(data_set_positive)
    result.loc["darnell_human_ViennaDuplex_features", "size"] = int(pos.shape[0])

    for method_dir in dir.iterdir():
        print(method_dir)
        for dataset_file_neg in method_dir.glob("*features*"):
            neg = read_csv(dataset_file_neg)

            result.loc[dataset_file_neg.name, "size"] = int(neg.shape[0])

    result = result.astype('int')
    print(result)
    summury = MERGE_DATA / "file size.csv"
    to_csv(result, summury)

split_train_test()
size()

