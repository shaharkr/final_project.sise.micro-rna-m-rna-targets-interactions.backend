from functools import partial
from pathlib import Path
from typing import Dict
import pandas as pd
import yaml
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.model_selection import PredefinedSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
import sklearn
from sklearn.svm import SVC
from xgboost import XGBClassifier
import pickle
import Classifier.FeatureReader as FeatureReader
from Classifier.FeatureReader import get_reader
from Classifier.ClfLogger import logger
from consts.global_consts import  ROOT_PATH, DATA_PATH_INTERACTIONS, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS

class ClassifierWithGridSearch(object):
    def __init__(self, dataset_file, result_dir):
        self.dataset_file = dataset_file
        self.dataset_name = self.extract_dataset_name()
        print(f"Handling dataset : {self.dataset_name}")
        self.load_dataset()
        self.result_dir = Path(result_dir)
        self.result_dir.mkdir(exist_ok=True, parents=True)
        self.create_clf_dict()

    def extract_dataset_name(self):
        return str(self.dataset_file.stem).split("_test")[0]

    def create_clf_dict(self):
        self.clf_dict = {
            "xgbs": XGBClassifier(),
            "xgbs_no_encoding": XGBClassifier(),
        }

    # this function response on load the dataset
    def load_dataset(self):
        directory = self.dataset_file.parent
        feature_reader = get_reader()
        X, y = feature_reader.file_reader(directory / f"{self.dataset_name}.csv")
        self.X = X
        # self.X.drop(columns=['MRNA_Target_A_comp' ,'MRNA_Target_C_comp','MRNA_Target_G_comp',
        # 'MRNA_Target_U_comp','MRNA_Target_AA_comp','MRNA_Target_AC_comp','MRNA_Target_AG_comp','MRNA_Target_AU_comp',
        # 'MRNA_Target_CA_comp','MRNA_Target_CC_comp','MRNA_Target_CG_comp','MRNA_Target_CU_comp',
        # 'MRNA_Target_GA_comp','MRNA_Target_GC_comp','MRNA_Target_GG_comp','MRNA_Target_GU_comp',
        # 'MRNA_Target_UA_comp','MRNA_Target_UC_comp','MRNA_Target_UG_comp','MRNA_Target_UU_comp'], inplace=True)
        self.y = y

    # This function response on train model and then save this model
    def train_one_conf_kfold(self, clf_name, conf, scoring="accuracy"):
        return_dict={}
        # open the path of result file
        output_file = self.result_dir / f"{self.dataset_name}_{clf_name}.csv"
        if output_file.is_file():
            print(f"output file: {output_file} exits. skip.")
            return

        # creat the specific clf and load the parameters of the clf according to the ymal file.
        clf = self.clf_dict[clf_name]
        print(clf)
        parameters = conf['parameters']
        grid_obj = GridSearchCV(estimator=clf,param_grid=parameters, scoring=scoring, cv=5, n_jobs=-1, verbose=3)
        # grid_obj = GridSearchCV(estimator=clf,param_grid= parameters, scoring=scoring, cv=1, n_jobs=-1, verbose=3)

        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=44)

        for train_index, test_index in skf.split(self.X, self.y):
            x0, x1 = self.X.iloc[train_index], self.X.iloc[test_index]
            y0, y1 = self.y[train_index], self.y[test_index]
            grid_obj.fit(x0, y0, eval_set=[(x1, y1)],
                    eval_metric='logloss', verbose=False, early_stopping_rounds=10)
            prval = grid_obj.predict_proba(x1)[:, 1]
            return_dict[self.dataset_name] = roc_auc_score(y1, prval)
        # this step find to optimize prams
        grid_obj = GridSearchCV(clf, parameters, scoring=scoring, cv=5, n_jobs=-1, verbose=3)
        grid_obj.fit(self.X, self.y)
        print('\n Best estimator:')
        print(grid_obj.best_estimator_)
        print(grid_obj.best_score_ * 2 - 1)
        # save the best classifier
        best_clf = grid_obj.best_estimator_
        model_file = self.result_dir / f"{self.dataset_name}_{clf_name}.model"

        try:
            with model_file.open("wb") as pfile:
                pickle.dump(best_clf, pfile)
        except Exception:
            pass
        results = pd.DataFrame(grid_obj.cv_results_)
        results.to_csv(output_file, index=False)

    def fit_best_clf(self, clf_name):
        clf = self.clf_dict[clf_name]
        print(clf)
        fit_params = {}
        if clf_name =="xgbs":
            fit_params = {"eval_set": [(self.X, self.y)],
                          "early_stopping_rounds": 50}

        clf.fit(self.X, self.y, **fit_params)
        clf.save_model("model_sklearn.json")

        return clf

    # this function response on train model and then save this model
    def train_one_conf(self, clf_name, conf, scoring="accuracy"):

        output_file = self.result_dir / f"{self.dataset_name}_{clf_name}.csv"
        # if output_file.is_file():
        #     print(f"output file: {output_file} exits. skip.")
        #     return

        # creat the specific clf and load the parameters of the clf according to the ymal file.
        clf = self.clf_dict[clf_name]
        print(clf)
        parameters = conf['parameters']

        # this step find to optimize prams
        grid_obj = GridSearchCV(clf, parameters, scoring=scoring, cv=5, n_jobs=-1, verbose=3)
        grid_obj.fit(self.X, self.y)

        print('\n Best estimator:')
        print(grid_obj.best_estimator_)
        print(grid_obj.best_score_ * 2 - 1)
        # save the best classifier
        best_clf = grid_obj.best_estimator_
        model_file = self.result_dir / f"{self.dataset_name}_{clf_name}.model"

        try:
            with model_file.open("wb") as pfile:
                pickle.dump(best_clf, pfile)
        except Exception:
            pass
        results = pd.DataFrame(grid_obj.cv_results_)
        results.to_csv(output_file, index=False)

        def fit_best_clf(self, clf_name):
            clf = self.clf_dict[clf_name]
            print(clf)
            fit_params = {}
            if clf_name == "xgbs":
                fit_params = {"eval_set": [(self.X, self.y)],
                              "early_stopping_rounds": 50}

            clf.fit(self.X, self.y, **fit_params)
            clf.save_model("model_sklearn.json")

            return clf

    def fit(self, yaml_path):
        with open(yaml_path, 'r') as stream:
            training_config = yaml.safe_load(stream)

        for clf_name, conf in training_config.items():
            if conf["run"]:
                self.train_one_conf(clf_name, conf, scoring="accuracy")
                #clf = self.fit_best_clf(clf_name)
                #print(clf)


def worker(dataset_file, results_dir, yaml_file):
    clf_grid_search = ClassifierWithGridSearch(dataset_file=dataset_file, result_dir=results_dir)
    clf_grid_search.fit(yaml_file)
    return


def self_fit(feature_mode, yaml_file, first_self, last_self, name_method, dir_method,number_iteration):
    logger.info("starting self_fit")
    logger.info(f"params: {[feature_mode, yaml_file, first_self, last_self]}")

    FeatureReader.reader_selection_parameter = feature_mode
    csv_dir = DATA_PATH_INTERACTIONS / "train" / name_method / number_iteration
    files = list(csv_dir.glob('**/*.csv'))
    for f in files:
        # if "non_overlapping_sites_darnell_human_ViennaDuplex_negative_features_train_underSampling_method_0" not in f.name:
        #     continue
        results_dir = ROOT_PATH / "Results/models" / dir_method / number_iteration
        logger.info(f"results_dir = {results_dir}")
        logger.info(f"start dataset = {f}")
        worker(f, results_dir=results_dir, yaml_file=yaml_file)
        logger.info(f"finish dataset = {f}")
    logger.info("finish self_fit")


def build_classifiers(number_iteration):
    yaml_file = "/sise/home/efrco/efrco-master/Classifier/yaml/xgbs_params_small.yml"
    # yaml_file = "/sise/home/efrco/efrco-master/Classifier/yaml/xgbs_params.yml"

    FeatureReader.reader_selection_parameter = "without_hot_encoding"
    number_iteration = str(number_iteration)
    self_fit("without_hot_encoding", yaml_file, 1, 2, name_method="underSampling", dir_method="models_underSampling", number_iteration=number_iteration)
    print("END main_primary")

# build_classifiers()

# build_classifiers(number_iteration=0)





