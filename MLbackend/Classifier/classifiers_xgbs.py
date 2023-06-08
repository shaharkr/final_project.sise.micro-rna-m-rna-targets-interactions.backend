from functools import partial
from pathlib import Path
from typing import Dict

import click
import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.model_selection import PredefinedSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import RandomizedSearchCV,GridSearchCV
from multiprocessing import Process, Pool
import sklearn

from sklearn.svm import SVC
from xgboost import XGBClassifier
import pickle
import FeatureReader
from FeatureReader import get_reader, reader_dict
from ClfLogger import logger




class ClassifierWithGridSearch (object):
    def __init__(self, dataset_file, result_dir):
        self.dataset_file = dataset_file
        self.dataset_name = self.extract_dataset_name()
        print (f"Handling dataset : {self.dataset_name}")
        self.load_dataset()
        self.result_dir = Path(result_dir)
        self.result_dir.mkdir(exist_ok=True, parents=True)
        self.create_clf_dict()

    def extract_dataset_name (self):
        return str(self.dataset_file.stem).split("_test")[0]

    def create_clf_dict (self):
        self.clf_dict = {
            'rf': RandomForestClassifier(),
            'SVM': SVC(),
            'logit': LogisticRegression(),
            'SGD': SGDClassifier(),
            "KNeighborsClassifier": KNeighborsClassifier(),
            "xgbs": XGBClassifier(),
            "xgbs_no_encoding" : XGBClassifier(),
        }


    def load_dataset(self):
        directory = self.dataset_file.parent
        feature_reader = get_reader()
        X, y = feature_reader.file_reader(directory / f"{self.dataset_name}_train.csv")
        self.X = X
        self.y = y


    def train_one_conf (self, clf_name, conf, scoring="accuracy"):
        output_file = self.result_dir / f"{self.dataset_name}_{clf_name}.csv"
        if output_file.is_file():
            print (f"output file: {output_file} exits. skip.")
            return

        clf = self.clf_dict[clf_name]
        print (clf)
        parameters = conf['parameters']

        grid_obj = GridSearchCV(clf, parameters, scoring=scoring, cv=4, n_jobs=-1, verbose=3)
        grid_obj.fit(self.X, self.y)

        print('\n Best estimator:')
        print(grid_obj.best_estimator_)
        print(grid_obj.best_score_ * 2 - 1)
        #save the best classifier
        best_clf = grid_obj.best_estimator_
        model_file = self.result_dir / f"{self.dataset_name}_{clf_name}.model"

        try:
            with model_file.open("wb") as pfile:
                pickle.dump(best_clf, pfile)
        except Exception:
            pass


        results = pd.DataFrame(grid_obj.cv_results_)
        results.to_csv(output_file, index=False)

    def fit_best_clf(self, clf_name, parameters: Dict):
        clf = self.clf_dict[clf_name]
        print (clf)
        fit_params = {}
        if clf_name=="xgbs":
            fit_params = {"eval_set": [(self.X_val, self.y_val)],
                          "early_stopping_rounds": 50}

        clf.fit(self.X, self.y, **fit_params)


        return clf

    def fit (self, yaml_path):
        with open(yaml_path, 'r') as stream:
            training_config = yaml.safe_load(stream)

        for clf_name, conf in training_config.items():
            if conf["run"]:
                self.train_one_conf(clf_name, conf, scoring="accuracy")


def worker(dataset_file, results_dir, yaml_file):
    clf_grid_search = ClassifierWithGridSearch(dataset_file=dataset_file, result_dir=results_dir)
    clf_grid_search.fit(yaml_file)
    return




def self_fit(feature_mode, yaml_file, first_self, last_self):
    logger.info("starting self_fit")
    logger.info(f"params: {[feature_mode, yaml_file, first_self, last_self]}")

    FeatureReader.reader_selection_parameter = feature_mode
    csv_dir = Path("Features/CSV")
    for i in range(first_self, last_self):
        train_test_dir = csv_dir / f"train_test{i}"
        results_dir = Path("Results") / f"self{i}"
        logger.info(f"train_test_dir = {train_test_dir}")
        logger.info(f"results_dir = {results_dir}")
        for dataset_file in train_test_dir.glob("*_test*"):
            logger.info(f"start dataset = {dataset_file}")
            worker(dataset_file, results_dir=results_dir, yaml_file=yaml_file)
            logger.info(f"finish dataset = {dataset_file}")
    logger.info("finish self_fit")





def different_fit(yaml_file, first_self, last_self):
    csv_dir = Path("Features/CSV")
    for i in range(first_self, last_self):
        train_test_dir = csv_dir / f"train_test{i}"
        results_dir = Path("Results") / f"self{i}" / "xgbs_different"
        results_dir.mkdir(exist_ok=True)
        for dataset_file in train_test_dir.glob("*test*"):
            worker(dataset_file, results_dir, yaml_file, True)




