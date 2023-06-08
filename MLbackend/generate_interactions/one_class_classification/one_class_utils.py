from consts.global_consts import  ROOT_PATH, DATA_PATH_INTERACTIONS, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from utilsfile import read_csv, to_csv
import Classifier.FeatureReader as FeatureReader
from Classifier.FeatureReader import get_reader
import pandas as pd
from collections import Counter


def split_train_one_class(method_split_source, random_state,  number_split):
    csv_dir = DATA_PATH_INTERACTIONS / "train" / method_split_source
    files = list(csv_dir.glob('**/*.csv'))
    dir_results = DATA_PATH_INTERACTIONS / "train" / 'one_class_svm' / str(number_split)

    for f in files:
        cols = []
        X_test = read_csv(f)
        # pos = X_test[X_test['Label'] == 1]
        # neg = X_test[X_test['Label'] == 0]
        # # print(pos.shape)
        # print("before:" , neg.shape)
        # neg = neg[neg.isna().any(axis=1)]
        # print("after:" , neg.shape)

        # cols = neg.columns[neg.isna().any()].tolist()
        # print(cols)
        # path = "/sise/home/efrco/efrco-master/mockmrnanull.csv"
        # to_csv(neg, path)


        df_new = X_test.drop(X_test[X_test['Label'] == 0].sample(frac=1, random_state=random_state).index)
        dataset = str(f.stem).split('/')[-1].split('_features')[0] + "_train_one_class.csv"
        out_results = dir_results / dataset

        # df_new = df_new.dropna()
        # pos = df_new[df_new['Label'] == 1]
        # neg = df_new[df_new['Label'] == 0]
        # print(pos.shape)
        # print(neg.shape)
        print("#################################################")
        # df_new.drop(columns=["Label"], inplace=True)
        to_csv(df_new, out_results)

def split_test_one_class(method_split_source, random_state,  number_split):
    csv_dir = DATA_PATH_INTERACTIONS / "test" / method_split_source
    files = list(csv_dir.glob('**/*.csv'))
    dir_results = DATA_PATH_INTERACTIONS / "test" / 'one_class_svm' / str(number_split)
    FeatureReader.reader_selection_parameter = "without_hot_encoding"

    for f in files:
        X_test = read_csv(f)
        # pos = X_test[X_test['Label'] == 1]
        # neg = X_test[X_test['Label'] == 0]
        # # print(pos.shape)
        # print(neg.shape)
        # neg = neg[neg.isna().any(axis=1)]
        # print(neg.shape)

        # cols = neg.columns[neg.isna().any()].tolist()
        # path = "/sise/home/efrco/efrco-master/enull.csv"
        # to_csv(C1, path)

        print("###########################")
        df_new = X_test.drop(X_test[X_test['Label'] == 0].sample(frac=.8, random_state=random_state).index)
        dataset = str(f.stem).split('/')[-1].split('_features')[0] + "_test_one_class.csv"
        out_results = dir_results / dataset
        # df_new.drop(columns=["Label"], inplace=True)
        # df_new = df_new.dropna()
        # pos = df_new[df_new['Label'] == 1]
        # neg = df_new[df_new['Label'] == 0]
        # print(pos.shape)
        # print(neg.shape)
        to_csv(df_new, out_results)

def split_train_test_one_class(method_split_source, random_state, number_split):
    split_train_one_class(method_split_source=method_split_source, random_state=random_state, number_split=number_split)
    split_test_one_class(method_split_source= method_split_source, random_state=random_state, number_split=number_split)





# def split_test_one_class(method_split_source, random_state):
#     csv_dir = DATA_PATH_INTERACTIONS / "test" / method_split_source
#     files = list(csv_dir.glob('**/*.csv'))
#     dir_results = DATA_PATH_INTERACTIONS / "test" / 'one_class_svm'
#     FeatureReader.reader_selection_parameter = "without_hot_encoding"
#     feature_reader = get_reader()
#
#     for f in files:
#         X_test = read_csv(f)
#         y_test = X_test.Label.ravel()
#         X_test.drop(columns=["Label"], inplace =True)
#         df_new = pd.concat([pd.DataFrame(y_test), X_test], axis=1)
#         df_new.rename({0: "Label"}, axis='columns', inplace=True)
#         df_new = df_new.dropna()
#         df_new = df_new.drop(df_new[df_new['Label'] == 0].sample(frac=.8, random_state=random_state).index)
#         dataset = str(f.stem).split('/')[-1].split('_features')[0] + "_test_one_class.csv"
#         out_results = dir_results / dataset

        # to_csv(df_new, out_results)

# split_train_one_class(method_split_source='underSampling', random_state=1, number_split=1)
# split_test_one_class(method_split_source='underSampling', random_state=1, number_split=1)
# split_train_one_class(method_split_source='underSampling', random_state=1)


def test_size():

    print("###############TEST###################")
    dir = DATA_PATH_INTERACTIONS / "test/one_class_svm"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*test*"):
        dataset = str(dataset_file.stem).split("_test")[0]
        print('#####################################################################')
        print(dataset)
        for suffix in ["test_one_class"]:
            dir = DATA_PATH_INTERACTIONS / "test/one_class_svm"
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
    dir = DATA_PATH_INTERACTIONS / "train/one_class_svm"
    result = pd.DataFrame()
    for dataset_file in dir.glob("*train*"):
        dataset = str(dataset_file.stem).split("_train")[0]
        print('#####################################################################')
        print(dataset)
        for suffix in ["train_one_class"]:
            dir = DATA_PATH_INTERACTIONS / "train/one_class_svm"
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