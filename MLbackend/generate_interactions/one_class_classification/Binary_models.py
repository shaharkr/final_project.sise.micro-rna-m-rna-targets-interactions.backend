from functools import partial
from pathlib import Path
from typing import Dict
import pandas as pd
from sklearn.datasets import make_classification
# Data processing
import pandas as pd
import yaml
import numpy as np
from collections import Counter
# Visualization
import matplotlib.pyplot as plt
# Model and performance
from sklearn.svm import OneClassSVM
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
import sklearn
from sklearn.svm import SVC
from xgboost import XGBClassifier
import pickle
from utilsfile import read_csv, to_csv
import Classifier.FeatureReader as FeatureReader
from Classifier.FeatureReader import get_reader
from sklearn.ensemble import RandomForestClassifier
from Classifier.ClfLogger import logger
from consts.global_consts import  ROOT_PATH, DATA_PATH_INTERACTIONS, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from sklearn.ensemble import IsolationForest

class NoModelFound(Exception):
    pass

class Binary_comper_model(object):
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
            "SVM":  SVC(),
            "rf": RandomForestClassifier()

        }

    # this function response on load the dataset
    def load_dataset(self):
        directory = self.dataset_file.parent
        feature_reader = get_reader()
        X, y = feature_reader.file_reader(directory / f"{self.dataset_name}.csv")
        self.X = X
        self.y= y


    # this function response on train model and then save this model
    def train_one_conf(self, clf_name, conf, scoring="accuracy"):

        output_file = self.result_dir / f"{self.dataset_name}_{clf_name}.csv"

        # creat the specific clf and load the parameters of the clf according to the ymal file.
        clf = self.clf_dict[clf_name]
        print(clf)
        parameters = conf['parameters']

        # this step find to optimize prams
        grid_obj = GridSearchCV(clf, parameters, scoring=scoring, cv=4, n_jobs=-1, verbose=3)
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

    def fit(self, yaml_path):
        with open(yaml_path, 'r') as stream:
            training_config = yaml.safe_load(stream)

        for clf_name, conf in training_config.items():
            if conf["run"]:
                self.train_one_conf(clf_name, conf, scoring="accuracy")


def worker(dataset_file, results_dir, yaml_file):
    clf_grid_search = Binary_comper_model(dataset_file=dataset_file, result_dir=results_dir)
    clf_grid_search.fit(yaml_file)
    return


def self_fit(feature_mode, yaml_file, first_self, last_self, name_method, dir_method):
    logger.info("starting self_fit")
    logger.info(f"params: {[feature_mode, yaml_file, first_self, last_self]}")

    FeatureReader.reader_selection_parameter = feature_mode
    csv_dir = DATA_PATH_INTERACTIONS / "train" / name_method
    files = list(csv_dir.glob('**/*.csv'))
    for f in files:
        results_dir = ROOT_PATH / "Results/models" / dir_method
        logger.info(f"results_dir = {results_dir}")
        logger.info(f"start dataset = {f}")
        worker(f, results_dir=results_dir, yaml_file=yaml_file)
        logger.info(f"finish dataset = {f}")
    logger.info("finish self_fit")


def build_classifiers_binary_comper():
    yaml_file = "/sise/home/efrco/efrco-master/Classifier/yaml/Binary_comper_model.yml"
    FeatureReader.reader_selection_parameter = "without_hot_encoding"
    self_fit("without_hot_encoding", yaml_file, 1, 2, name_method="underSampling", dir_method="models_binary_comper")
    print("END main_primary")

build_classifiers_binary_comper()
