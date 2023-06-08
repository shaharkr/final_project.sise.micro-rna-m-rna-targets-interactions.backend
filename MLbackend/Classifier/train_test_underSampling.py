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
from utils.utilsfile import read_csv, to_csv
# from utilsfile import read_csv, to_csv



def imbalanced_partation(pos, pos_train, pos_test, input_file_negative, output_dir, test_size, random_state, number_split):

    neg = read_csv(input_file_negative)
    neg.insert(0, "Label", 0)
    number_split = str(number_split)
    # Both dataset must have the same columns
    col = [c for c in pos.columns if c in neg.columns]
    pos_test = pos_test[col]
    pos_train = pos_train[col]
    pos = pos[col]
    neg = neg[col]
    if neg.shape[0] > 100000:
        neg = neg.sample(n=100000, random_state=random_state)

    neg.insert(0, "Label1", 0)
    pos.insert(0, "Label1", 1)

    merge = pd.concat([pos, neg])
    y = merge.Label.ravel()
    X = merge.drop(columns=['Label'], axis=1)
    print(Counter(y))
    # define undersample strategy
    undersample = RandomUnderSampler(sampling_strategy='majority', random_state=random_state)

    # fit and apply the transform
    X_over, y_over = undersample.fit_resample(X, y)
    print(Counter(y_over))
    X_over = X_over.rename(columns={"Label1": "Label"})

    x_over_negative = X_over[X_over['Label'] == 0]
    train, test = train_test_split(x_over_negative, test_size=0.2,  random_state=random_state)

    output = {
        "train": pos_train.append(train, ignore_index=True),
        "test": pos_test.append(test, ignore_index=True)
    }

    # save to csv
    output_dir = DATA_PATH_INTERACTIONS
    dataset = str(input_file_negative.stem).split("_features_")[0]
    for key, d in output.items():
        out_file = f"{dataset}_{key}"
        out_file = out_file + f"_underSampling_method_{number_split}.csv"
        d = d.reindex(np.random.RandomState(seed=random_state).permutation(d.index))
        output_dir_train_test = output_dir / key / "underSampling" / number_split
        to_csv(d, output_dir_train_test / out_file)
        del(d)



# function that response to split with random
def split_train_test(dataset_positive_name, random_state, number_split):

    dir = NEGATIVE_DATA_PATH
    train_test_dir = DATA_PATH_INTERACTIONS

    # Here we choose one dataset from all the data set of positive
    data_set_positive = MERGE_DATA /"positive_interactions_new/featuers_step"/ f"{dataset_positive_name}.csv"

    # first we dived the positive interactions, the partition is const for all the train and test set.
    print(f"working on file: pos={data_set_positive}")
    pos = read_csv(data_set_positive)
    pos.insert(0, "Label", 1)
    print("pos", data_set_positive, pos.shape)

    pos_train, pos_test = train_test_split(pos, test_size=0.2, random_state=random_state)
    for method_dir in dir.iterdir():
        list_method = ['mockMirna', 'non_overlapping_sites']
        for dataset_file_neg in method_dir.glob("*features*"):
            # Only take the correct mock file
            if any(method in method_dir.stem for method in list_method):
                correct_dataset = dataset_positive_name.split("_features")[0] + "_negative_features"
                list_split = dataset_file_neg.stem.split("/")
                list_split = [str(s) for s in list_split]
                if correct_dataset not in list_split[0]:
                    continue
            imbalanced_partation(pos, pos_train, pos_test, dataset_file_neg, train_test_dir, 0.2, random_state, number_split)

# split_train_test(dataset_positive_name="darnell_human_ViennaDuplex_features", random_state=19, number_split=0)

def test_size():

    print("###############TEST###################")
    dir = DATA_PATH_INTERACTIONS / "test/underSampling"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*test*"):
        dataset = str(dataset_file.stem).split("_test")[0]
        print('#####################################################################')
        print(dataset)
        for suffix in ["test_underSampling_method"]:
            dir = DATA_PATH_INTERACTIONS / "test/underSampling"
            d = read_csv(dir / f"{dataset}_{suffix}.csv")
            suffix = suffix.replace("test_", "")
            result.loc[dataset, suffix] = int(d.shape[0])
            print(Counter(d['Label']))
    result.sort_index(inplace=True)
    result = result.astype('int')
    print(result)
    print(60 * "*")
    # print(result.to_latex())


def train_size():

    print("###############Train###################")
    dir = DATA_PATH_INTERACTIONS / "train/underSampling"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*train*"):
        dataset = str(dataset_file.stem).split("_train")[0]
        print('#####################################################################')
        print(dataset)
        for suffix in ["train_underSampling_method"]:
            dir = DATA_PATH_INTERACTIONS / "train/underSampling"
        d = read_csv(dir / f"{dataset}_{suffix}.csv")
        suffix = suffix.replace("train_", "")
        result.loc[dataset, suffix] = int(d.shape[0])
        print(Counter(d['Label']))

    result.sort_index(inplace=True)
    result = result.astype('int')
    print(result)
    print(60 * "*")
    # print(result.to_latex())

#
# test_size()
# train_size()
# split_train_test("darnell_human_ViennaDuplex_features", random_state=0, number_split=0)
