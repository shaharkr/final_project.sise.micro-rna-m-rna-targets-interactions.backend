from pathlib import Path
import sys
sys.path.append('/sise/home/efrco/efrco-master/Classifier/')
from Classifier.train_test_underSampling import split_train_test as split_train_test_underSampling
from Classifier.ClassifierWithGridSearch import build_classifiers as build_classifiers_grid_search
from Classifier.result_test import different_results_summary
import pandas as pd
from utils.utilsfile import read_csv, to_csv
from consts.global_consts import ROOT_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from generate_interactions.one_class_classification.one_class import build_classifiers_svm
from generate_interactions.one_class_classification.Isolation_Forest import build_classifiers_isolation_forest
from generate_interactions.one_class_classification.one_class_utils import split_test_one_class, split_train_one_class

class NoModelFound(Exception):
    pass


class CrossValidation(object):

    def __init__(self, dataset_file_positive, result_dir, number_iterations):
        self.number_iterations = number_iterations
        self.dataset_file_positive = dataset_file_positive
        self.result_dir = ROOT_PATH / Path("Results") / result_dir
        self.result_dir.mkdir(exist_ok=True, parents=True)
        # self.measurement_dict_xgboost = self.screate_measurement_dict()
        # self.measurement_dict_svm = self.create_measurement_dict()
        self.measurement_dict ={'svm': self.create_measurement_dict(),
                                'xgbs': self.create_measurement_dict(),
                                'isolation_forest': self.create_measurement_dict()}
        # self.measurement_dict = {
        #                          'isolation_forest': self.create_measurement_dict()}


    def create_measurement_dict(self):
        columns_list = [f"exper_{i}" for i in range(self.number_iterations)]
        measurement_dict = {
            "AUC": pd.DataFrame(columns=columns_list, dtype=object),
            "ACC": pd.DataFrame(columns=columns_list, dtype=object),
            "TPR": pd.DataFrame(columns=columns_list, dtype=object),
            "TNR": pd.DataFrame(columns=columns_list, dtype=object),
            "PPV": pd.DataFrame(columns=columns_list, dtype=object),
            "NPV": pd.DataFrame(columns=columns_list, dtype=object),
            "FPR": pd.DataFrame(columns=columns_list, dtype=object),
            "FNR": pd.DataFrame(columns=columns_list, dtype=object),
            "FDR": pd.DataFrame(columns=columns_list, dtype=object),
            "F1": pd.DataFrame(columns=columns_list, dtype=object),

        }
        return measurement_dict

    def get_measurement_dict(self, name_dict):
        return {k: round(v, 3) for k, v in self.measurement_dict[name_dict].items()}



    def run_xgboost(self, number_iteration):
        build_classifiers_grid_search()
        results_dir = ROOT_PATH / Path("Results")
        different_results_summary(method_split="underSampling", model_dir="models_underSampling", number_iteration=number_iteration, name_classifier='xgbs')
        ms_table = read_csv(results_dir / 'results_iterations'/ "xgbs" / f"measurement_summary_{number_iteration}.csv")
        for measurement in self.measurement_dict['xgbs'].keys():
            col = ms_table[measurement].apply(lambda t: round(t, 3))
            self.measurement_dict['xgbs'][measurement][f"exper_{number_iteration}"] = col

        # save file of result for each measuerment

    def run_one_class_svm(self, number_iteration):
        build_classifiers_svm()
        different_results_summary(method_split="one_class_svm", model_dir="models_one_class_svm",
                                  number_iteration=number_iteration, name_classifier='svm')
        results_dir = ROOT_PATH / Path("Results")
        ms_table = read_csv(results_dir / 'results_iterations' / "svm" / f"measurement_summary_{number_iteration}.csv")
        for measurement in self.measurement_dict['svm'].keys():
            col = ms_table[measurement].apply(lambda t: round(t, 3))
            self.measurement_dict['svm'][measurement][f"exper_{number_iteration}"] = col

    def run_isolation_forest(self, number_iteration):
        build_classifiers_isolation_forest()
        different_results_summary(method_split="one_class_svm", model_dir="models_isolation_forest",
                                  number_iteration=number_iteration, name_classifier='isolation_forest')
        results_dir = ROOT_PATH / Path("Results")
        ms_table = read_csv(results_dir / 'results_iterations' / "isolation_forest" / f"measurement_summary_{number_iteration}.csv")
        for measurement in self.measurement_dict['isolation_forest'].keys():
            col = ms_table[measurement].apply(lambda t: round(t, 3))
            self.measurement_dict['isolation_forest'][measurement][f"exper_{number_iteration}"] = col


    def write_results(self):

        # save file of result for each measuerment
        for classifier in self.measurement_dict.keys():
            for measurement in self.get_measurement_dict(classifier).keys():
                out_dir = self.result_dir/classifier/f"{measurement}_summary.csv"
                to_csv(self.get_measurement_dict(classifier)[measurement], out_dir)

    def run_experiment(self, method_split, model_dir):

        for i in range(self.number_iterations):
            split_train_test_underSampling(dataset_positive_name=self.dataset_file_positive, random_state=i*19)
            self.run_xgboost(number_iteration=i)
            split_train_one_class(method_split_source='underSampling', random_state=i*19)
            split_test_one_class(method_split_source='underSampling', random_state=i*19)
            self.run_one_class_svm(number_iteration=i)
            self.run_isolation_forest(number_iteration=i)
        self.write_results()

    def clean_name(self, name):

        name_clean = name.replace("model:", "").replace("human", "").replace("darnell_human", "").replace("test:",
                                                                                                          "").replace(
            "ViennaDuplex", "").replace("_darnell_", "").replace("__", "")

        name_clean = name_clean.replace("_nucleotides", "_mono", )
        name_clean = name_clean.replace("denucleotides", "_de")
        name_clean = name_clean.replace("method1", "mrna")
        name_clean = name_clean.replace("method2", "site")
        name_clean = name_clean.replace("method3", "mockMirna")

        if name_clean == "mockMirnadarnell":
            name_clean = "mockMirna"
        if name_clean == "mockMrnadarnell":
            name_clean = "mockMrna"
        if name_clean == "nonoverlappingsitesdarnell":
            name_clean = "Non site"
        return name_clean

    def summary_matrix(self, measurement_name):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir/classifier / "final_mean.csv"
            df_mean = read_csv(out_dir_mean)
            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'mono' in train_name or 'mono' in test_name:
                    continue
                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]
            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")

            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/") / "heatmap.png"
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()


    def summary_matrix_mock_mrna(self, measurement_name):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir/classifier / "final_mean.csv"
            df_mean = read_csv(out_dir_mean)
            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'mockMrna' not in train_name or 'mockMrna' not in test_name:
                    continue

                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]
            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")

            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/") / "heatmap_mockmrna.png"
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_result(self):

        all_result_mean = pd.DataFrame()
        all_result_std = pd.DataFrame()

        for classifier in self.measurement_dict.keys():
            all_result_mean = pd.DataFrame()
            all_result_std = pd.DataFrame()
            dir_measurement = self.result_dir / classifier
            for measuerment_file in dir_measurement.glob("*summary*"):
                df = read_csv(measuerment_file)
                count = 0
                measuerment_name = measuerment_file.stem.split("_summary")[0]

                for index, row in df.iterrows():
                    row = df.iloc[count]
                    count = count + 1
                    col_mean = row.mean()
                    col_std = row.std()
                    all_result_mean.loc[measuerment_name, index] = round(col_mean, 3)
                    all_result_std.loc[measuerment_name, index] = round(col_std, 3)

            out_dir_mean = self.result_dir/classifier / f"final_mean.csv"
            to_csv(all_result_mean, out_dir_mean)

            out_dir_std = self.result_dir/classifier / f"final_std.csv"
            to_csv(all_result_std, out_dir_std)


############################### Runnning ############################################
cv = CrossValidation("darnell_human_ViennaDuplex_features", "measurments_cross_validation", number_iterations=1)
cv.run_experiment(method_split='underSampling', model_dir='models_underSampling')
cv.summary_result()
cv.summary_matrix('ACC')
cv.summary_matrix_mock_mrna('ACC')