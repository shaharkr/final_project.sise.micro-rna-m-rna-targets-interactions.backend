# from Classifier.train_test_split_stratify import split_train_test
from Classifier.train_test_underSampling import split_train_test

from Classifier.ClassifierWithGridSearch import build_classifiers
from Classifier.result_test import different_results_summary


def model_run():
    # split_train_test()
    # main_primary()
    # different_results_summary(method_split="validation", model_dir="models_validation")

    split_train_test()
    build_classifiers()
    different_results_summary(method_split="underSampling", model_dir="models_underSampling")


model_run()