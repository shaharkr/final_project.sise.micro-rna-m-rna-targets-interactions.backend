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
from generate_interactions.one_class_classification.one_class_utils import split_train_test_one_class
import os
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
        # self.measurement_dict ={'svm': self.create_measurement_dict(),
        #                         'xgbs': self.create_measurement_dict(),
        #                         'isolation_forest': self.create_measurement_dict()}
        self.measurement_dict = {
                                 'xgbs': self.create_measurement_dict()}


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


    def clean_directory(self):
        train = "/sise/home/efrco/efrco-master/data/train"
        test = "/sise/home/efrco/efrco-master/data/test"
        figuers = "/sise/home/efrco/efrco-master/Results/figuers/"
        results_iterations = "/sise/home/efrco/efrco-master/Results/results_iterations/"
        models_clean = "/sise/home/efrco/efrco-master/Results/results_iterations/"
        measurment_cross = "/sise/home/efrco/efrco-master/Results/measurments_cross_validation/"

        paths = [train, test, figuers, results_iterations, models_clean, measurment_cross]

        for p in paths:
            for dirpath, dirnames, filenames in os.walk(p):
                for filename in filenames:
                    file_path = os.path.join(dirpath, filename)
                    try:
                        os.remove(file_path)
                    except OSError:
                        print("Error while deleting file: ", file_path)



    def split_train_test_files(self):
        self.clean_directory()
        for i in range(self.number_iterations):
            split_train_test_underSampling(dataset_positive_name=self.dataset_file_positive, random_state=i * 19, number_split=i)
            split_train_test_one_class(method_split_source='underSampling', random_state=i * 19,  number_split=i)

    def run_xgboost(self, number_iteration):
        build_classifiers_grid_search(number_iteration)
        different_results_summary(method_split="underSampling", model_dir="models_underSampling", number_iteration=number_iteration, name_classifier='xgbs')


    def run_one_class_svm(self, number_iteration):
        build_classifiers_svm(number_iteration)
        different_results_summary(method_split="one_class_svm", model_dir="models_one_class_svm",
                                  number_iteration=number_iteration, name_classifier='svm')

    def run_isolation_forest(self, number_iteration):
        build_classifiers_isolation_forest(number_iteration)
        different_results_summary(method_split="one_class_svm", model_dir="models_isolation_forest",
                                  number_iteration=number_iteration, name_classifier='isolation_forest')

    def summary_results_do_dict(self):
        for classifier in self.measurement_dict.keys():
            for number_iteration in range(self.number_iterations):
                results_dir = ROOT_PATH / Path("Results")
                ms_table = read_csv(
                    results_dir / 'results_iterations' / classifier / f"measurement_summary_{number_iteration}.csv")
                for measurement in self.measurement_dict[classifier].keys():
                    col = ms_table[measurement].apply(lambda t: round(t, 3))
                    self.measurement_dict[classifier][measurement][f"exper_{number_iteration}"] = col

    def write_results(self):
        self.summary_results_do_dict()
        # save file of result for each measuerment
        for classifier in self.measurement_dict.keys():
            for measurement in self.get_measurement_dict(classifier).keys():
                out_dir = self.result_dir/classifier/f"{measurement}_summary.csv"
                to_csv(self.get_measurement_dict(classifier)[measurement], out_dir)
        self.summary_result()
        self.summary_matrix('ACC')
        self.summary_matrix_mock_mrna('ACC')

    def run_experiment_xgbs(self, start, to):
        for i in range(start, to):
            self.run_xgboost(number_iteration=i)

    def run_experiment_one_class_svm(self):
        for i in range(self.number_iterations):
           self.run_one_class_svm(number_iteration=i)

    def run_experiment_isolation_forest(self):
        for i in range(self.number_iterations):
            self.run_isolation_forest(number_iteration=i)

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
            sns.color_palette("Spectral", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")
            plt.xticks(rotation=30, ha='right')

            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/") / "heatmap.png"
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_barplot_inra_results(self, measurement_name):
        keys = ['method', 'measurement', 'value']
        # res_table = pd.DataFrame(columns=list(keys), dtype=object)
        dtypes = np.dtype(
            [
                ("method", str),
                ("measurement", str),
                ("value", float),

            ]
        )
        res_table = pd.DataFrame(np.empty(0, dtype=dtypes))

        for classifier in self.measurement_dict.keys():

            out_dir_mean = self.result_dir / classifier / "final_mean.csv"
            df_mean = read_csv(out_dir_mean)
            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'mono' in train_name or 'mono' in test_name:
                    continue
                if train_name != test_name:
                    continue
                res_table.loc['method'] = train_name
                res_table.loc['measurement'] = measurement_name
                res_table.loc['value'] = df_mean.loc[measurement_name, columnName]


            ax = sns.barplot(data=res_table, x="measurement_name", y="value", hue="method")

            sns.color_palette("Spectral", as_cmap=True)


            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/") / "heatmap.png"
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_matrix_mock_mrna(self, measurement_name, target_calculate):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir/classifier / "final_mean.csv"
            out_dir_std = self.result_dir / classifier / "final_std.csv"

            if target_calculate == "mean":
                df_mean = read_csv(out_dir_mean)
            else:
                df_mean = read_csv(out_dir_std)

            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'mockMrna' not in train_name or 'mockMrna' not in test_name:
                    continue

                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]
            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")

            name_file = "heatmap_mockMrna" + str(measurement_name) + "_" + str(target_calculate) + ".png"
            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/mockMrna") / name_file
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()


    def summary_matrix_tarBase(self, measurement_name, target_calculate="mean"):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir / classifier / "final_mean.csv"
            out_dir_std = self.result_dir / classifier / "final_std.csv"

            if target_calculate == "mean":
                df_mean = read_csv(out_dir_mean)
            else:
                df_mean = read_csv(out_dir_std)

            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'tarBase' not in train_name or 'tarBase' not in test_name:
                    continue

                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]

            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")

            name_file = "heatmap_tarBase" + str(measurement_name) + "_" + str(target_calculate) + ".png"
            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/tarBase") / name_file
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_matrix_non_overlapping_site(self, measurement_name, target_calculate="mean"):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir / classifier / "final_mean.csv"
            out_dir_std = self.result_dir / classifier / "final_std.csv"

            if target_calculate == "mean":
                df_mean = read_csv(out_dir_mean)
            else:
                df_mean = read_csv(out_dir_std)

            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'non_overlapping_sites' not in train_name or 'non_overlapping_sites' not in test_name:
                    continue

                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]

            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")

            name_file = "heatmap_non_overlapping_sites" + str(measurement_name) + "_" + str(target_calculate) + ".png"
            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/non_overlapping_sites") / name_file
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_matrix_clip(self, measurement_name, target_calculate="mean"):
        res_table = pd.DataFrame()
        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir / classifier / "final_mean.csv"
            out_dir_std = self.result_dir / classifier / "final_std.csv"

            if target_calculate == "mean":
                df_mean = read_csv(out_dir_mean)
            else:
                df_mean = read_csv(out_dir_std)

            for (columnName, columnData) in df_mean.iteritems():
                train_name = self.clean_name(columnName.split('/')[0])
                test_name = self.clean_name(columnName.split('/')[1])
                if 'clip_3' not in train_name or 'clip_3' not in test_name:
                    continue

                res_table.loc[test_name, train_name] = df_mean.loc[measurement_name, columnName]


            ax = sns.heatmap(res_table, annot=True)
            sns.color_palette("rocket", as_cmap=True)
            ax.set(xlabel="train", ylabel="test")


            name_file = "heatmap_clip" + str(measurement_name) + "_" + str(target_calculate) + ".png"
            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap/clip") / name_file
            plt.savefig(fname, format="PNG", bbox_inches='tight')
            plt.clf()

    def summary_matrix_Intra(self, measurement_name):
        res_table = pd.DataFrame()
        measurement_name_list = ['ACC', 'TPR', 'TNR', 'PPV']
        final_method = ['tarBase', "non_overlapping_sites", "mockMirna","mockMrna__de_site",
                        "clip_interaction_clip_3_", "non_overlapping_sites_clip_data"]

        for classifier in self.measurement_dict.keys():
            res_table = pd.DataFrame()
            out_dir_mean = self.result_dir / classifier / "final_mean.csv"

            df_mean = read_csv(out_dir_mean)
            for measurement_name in measurement_name_list:
                for (columnName, columnData) in df_mean.iteritems():
                    train_name = self.clean_name(columnName.split('/')[0])
                    test_name = self.clean_name(columnName.split('/')[1])
                    if train_name != test_name:
                        continue
                    if train_name not in final_method:
                        print(train_name)
                        continue
                    train_name = train_name.replace("__","_")

                    # res_table.loc[measurement_name, train_name] = df_mean.loc[measurement_name, columnName]
                    res_table.loc[train_name, measurement_name] = df_mean.loc[measurement_name, columnName]

            # ax = sns.heatmap(res_table, annot=True, linewidth=4.5,cmap=sns.cubehelix_palette(as_cmap=True))

            ax = sns.heatmap(res_table, annot=True, linewidth=(3.5,0))
            sns.color_palette("rocket", as_cmap=True)

            ax.set(xlabel="measurement", ylabel="datset")

            plt.xticks(rotation=30, ha='right')

            name_file = "heatmap_Intra.png"
            fname = ROOT_PATH / Path(f"Results/figuers/{classifier}/heatmap") / name_file
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

    def run_pipline(self):
        # self.split_train_test_files()
        # self.run_experiment_xgbs(start=0, to=5)
        # self.run_experiment_xgbs(start=5, to=9)
        # self.run_experiment_xgbs(start=9, to=12)
        # self.run_experiment_xgbs(start=12, to=15)

        # self.run_experiment_xgbs(start=2, to=6)
        # self.run_experiment_xgbs(start=6, to=9)
        # self.run_experiment_xgbs(start=9, to=12)
        # self.run_experiment_xgbs(start=12, to=15)
        # self.run_experiment_xgbs(start=15, to=18)
        # self.run_experiment_xgbs(start=18, to=20)

        # self.run_experiment_one_class_svm()
        # self.run_experiment_isolation_forest()
        self.write_results()


############################### Runnning ############################################
cv = CrossValidation("darnell_human_ViennaDuplex_features", "measurments_cross_validation", number_iterations=1)
# cv.run_pipline()
measurment_cross = "/sise/home/efrco/efrco-master/Results/figuers/xgbs/confuse_matrix/"
c = "/sise/home/efrco/efrco-master/Results/figuers/xgbs/feature_importance/"
paths = [measurment_cross, c]

for p in paths:
    for dirpath, dirnames, filenames in os.walk(p):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            try:
                os.remove(file_path)
            except OSError:
                print("Error while deleting file: ", file_path)
different_results_summary(method_split="underSampling", model_dir="models_underSampling",
                          number_iteration=1, name_classifier='xgbs')

# cv.summary_barplot_inra_results('ACC')


