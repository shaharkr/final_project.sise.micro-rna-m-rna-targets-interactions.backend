import ast
import shap
import pickle
import os
from collections import Counter
from itertools import combinations
from pathlib import Path
# from utilsClassifier import mean_std
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame
from seaborn import heatmap
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
from sklearn.svm import SVC
from xgboost import XGBClassifier
# import FeatureReader
# from FeatureReader import get_reader
# from dataset import Dataset
# from consts.global_consts import ROOT_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
# from dataset import Dataset
# from utils.utilsfile import read_csv, to_csv
from utilsfile import read_csv, to_csv
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import Classifier.FeatureReader as FeatureReader
from Classifier.FeatureReader import get_reader
from Classifier.ClfLogger import logger
from consts.global_consts import ROOT_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
import seaborn as sns
from sklearn.metrics import plot_confusion_matrix



class NoModelFound(Exception):
    pass

def model_confuse_matrix(x_test, y_test, model, model_name,s_org,d_org,name_classifiers):

    plot_confusion_matrix(model, x_test, y_test)
    plt.title(f"{s_org}_{d_org}_confuse_matrix_plot")
    fname = ROOT_PATH / Path(f"Results/figuers/{name_classifiers}/confuse_matrix/") / f"{s_org}_{d_org}.pdf"
    plt.savefig(fname, format="pdf", bbox_inches='tight')
    plt.show()
    plt.clf()



def measurement(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred).ravel()
    # cm = confusion_matrix(y_true, y_pred)

    # print(cm)
    # TP = cm[0][0]
    # FP = cm[0][1]
    # FN = cm[1][0]
    # TN = cm[1][1]
    TP = cm[3]
    FP = cm[1]
    FN = cm[2]
    TN = cm[0]
    print("TP:", TP)
    print("FP:", FP)
    print("FN:", FN)
    print("TN:", TN)
    try:
        auc = roc_auc_score(y_true, y_pred)
    except:
        auc = 0

    d = {
        # Sensitivity, hit rate, recall, or true positive rate
        "TPR": TP / (TP + FN),
        # Specificity or true negative rate
        "TNR": TN / (TN + FP),
        # Precision or positive predictive value
        "PPV": TP / (TP + FP),
        # Negative predictive value
        "NPV": TN / (TN + FN),
        # Fall out or false positive rate
        "FPR": FP / (FP + TN),
        # False negative rate
        "FNR": FN / (TP + FN),
        # False discovery rate
        "FDR": FP / (TP + FP),
        # Overall accuracy
        "ACC": (TP + TN) / (TP + FP + FN + TN),
        # roc auc
        "AUC": auc,
        "F1": f1_score(y_true, y_pred)
    }

    return {k: round(v,3) for k, v in d.items()}


# this function response on load the clf from Result dir
def get_presaved_clf(results_dir: Path, dataset: str, method: str):
    clf_file = results_dir / f"{dataset}_{method}.model"
    if not clf_file.is_file():
        raise NoModelFound(f"No model found: {clf_file}")
    with clf_file.open("rb") as f:
        clf = pickle.load(f)
    return clf


def xgbs_feature_importance(clf: XGBClassifier, X_train: DataFrame):
    def series2bin(d: pd.DataFrame, bins: list):
        for l in bins:
            d[f"bins{l}"] = d["rank"].apply(lambda x: int(x/l))
        return d

    feature_importances = pd.DataFrame(clf.feature_importances_, index=X_train.columns, columns=['importance'])
    feature_importances.sort_values('importance', ascending=False, inplace=True)
    feature_importances["rank"] = range(feature_importances.shape[0])
    feature_importances = series2bin(feature_importances, [10, 20, 50])

    # feature_importance_plot(feature_importances)
    return feature_importances


def plot_feature_importances(clf: XGBClassifier, X_train, s_org, d_org, name_classifiers):
    feat_imp = pd.DataFrame({'importance': clf.feature_importances_})
    feat_imp['feature'] = X_train.columns
    feat_imp.sort_values(by='importance', ascending=False, inplace=True)
    feat_imp = feat_imp.iloc[:7]

    feat_imp.sort_values(by='importance', inplace=True)
    feat_imp = feat_imp.set_index('feature', drop=True)
    feat_imp.plot.barh(title="title")
    plt.xlabel('Feature Importance Score')
    s_org = s_org.split('ViennaDuplex')[0]
    d_org = d_org.split('ViennaDuplex')[0]

    # fname = ROOT_PATH / Path("Results/figuers/feature_importance/") / f"{s_org}_{d_org}_summary_plot.pdf"
    fname = ROOT_PATH / Path(f"Results/figuers/{name_classifiers}/feature_importance/") / f"{s_org}_{d_org}_summary_plot.pdf"
    plt.savefig(fname, format="pdf", bbox_inches='tight')
    plt.show()
    plt.clf()


def model_shap_plot(test, model, model_name,s_org,d_org,name_classifiers, dependence_feature=None):

    shap_values = shap.TreeExplainer(model).shap_values(test.values.astype('float'))
    shap.summary_plot(shap_values, test, show=False, max_display=10, feature_names=test.columns)
    plt.title(f"{s_org}_{d_org}_summary_plot")
    fname = ROOT_PATH / Path(f"Results/figuers/{name_classifiers}/shap/") / f"{s_org}_{d_org}_summary_plot.pdf"
    plt.savefig(fname, format="pdf", bbox_inches='tight')
    plt.show()
    plt.clf()
    # if dependence_feature is not None:
    #     shap.dependence_plot(dependence_feature, shap_values, data,show=False)
    #     plt.title(f"{model_name}_{s_org}_{d_org}_dependence_plot")
    #     plt.savefig(os.path.join(MODELS_FEATURE_DEPENDENCE, f"{model_name}_{s_org}_{d_org}_dependence_plot.png"),
    #                 bbox_inches='tight')
    #     plt.clf()


#
# def self_results_summary(method_split: str):
#     ms_table = pd.DataFrame(columns={"TPR", "TNR", "PPV", "NPV", "FPR", "FNR", "FDR","ACC"})
#     results_dir = ROOT_PATH / Path("Results")
#     test_dir = DATA_PATH_INTERACTIONS / "test" / method_split
#     train_dir = DATA_PATH_INTERACTIONS / "train"
#
#     methods = ['xgbs']
#     res_table: DataFrame = pd.DataFrame()
#     FeatureReader.reader_selection_parameter = "without_hot_encoding"
#     feature_reader = get_reader()
#
#     for f_test in test_dir.glob("*test*"):
#         f_stem = f_test.stem
#         dataset = f_stem.split(".csv")[0]
#         for method in methods:
#             print(f"test: {f_test}, method: {method}")
#             try:
#                 dataset = dataset.replace("test", "train")
#                 clf = get_presaved_clf(results_dir, dataset, method)
#                 X_test, y_test = feature_reader.file_reader(f_test)
#                 model_shap_plot(X_test, clf, method, dataset, f_test,name_classifiers,  dependence_feature=None)
#                 test_score = accuracy_score(y_test, clf.predict(X_test))
#                 res_table.loc[dataset, method] = round(test_score, 3)
#
#                 print(res_table)
#                 name_summary_file = "summary_" + dataset + ".csv"
#                 res_table.to_csv(results_dir / "summary" / name_summary_file)
#                 if method in ["xgbs_no_encoding", "xgbs"]:
#                     feature_importance = xgbs_feature_importance(clf, X_test)
#                     to_csv(feature_importance, results_dir / "feature_importance" / f"feature_importance_{dataset}.csv")
#                     print("save feature importance file")
#                     ms = measurement(y_test, clf.predict(X_test))
#                     ms_table.loc[dataset, f_test] = ms
#
#             except NoModelFound:
#                 pass
#
#     res_table.sort_index(inplace=True)
#     print(res_table)
#     print(res_table.to_latex())
#
#     res_table.to_csv(results_dir/"summary.csv")
#     to_csv(ms_table, results_dir / "xgbs_measurements" / "xgbs_measurements.csv")
#

def different_results_summary(method_split: str, model_dir: str, number_iteration: int, name_classifier: str):

    ms_table = None
    results_dir = ROOT_PATH / Path("Results")
    number_iteration = str(number_iteration)
    results_dir_models = ROOT_PATH / Path("Results/models") / model_dir / number_iteration
    test_dir = DATA_PATH_INTERACTIONS / "test" / method_split / number_iteration
    res_table: DataFrame = pd.DataFrame()
    FeatureReader.reader_selection_parameter = "without_hot_encoding"
    feature_reader = get_reader()

    clf_datasets = [f.stem.split("_"+ name_classifier)[0] for f in results_dir_models.glob("*.model")]
    method = name_classifier
    for clf_dataset in clf_datasets:
        for f_test in test_dir.glob("*test*"):
            f_stem = f_test.stem
            test_dataset = f_stem.split(".csv")[0]
            print(f"clf: {clf_dataset} test: {test_dataset}, method: {method}")
            try:
                    clf = get_presaved_clf(results_dir_models, clf_dataset, method)
                    X_test, y_test = feature_reader.file_reader(test_dir/f"{test_dataset}.csv")

                    # score predict
                    test_score = accuracy_score(y_test, clf.predict(X_test))
                    res_table.loc[clf_dataset, test_dataset] = round(test_score, 3)

                    # feature_importance = xgbs_feature_importance(clf, X_test)
                    # name_method = "model:" + clf_dataset.split("negative")[0] + "_" + "test:" + test_dataset.split("negative")[0]
                    # to_csv(feature_importance, results_dir / "feature_importance" / f"feature_importance_{name_method}.csv")

                    # save measures
                    if name_classifier == 'svm' or name_classifier == 'isolation_forest':
                        # Get the scores for the testing dataset
                        score = clf.score_samples(X_test)
                        # Check the score for 2% of outliers
                        score_threshold = np.percentile(score, 2)
                        print(f'The customized score threshold for 2% of outliers is {score_threshold:.2f}')
                        # Check the model performance at 2% threshold
                        # prediction_1 = [0 if i < score_threshold else 1 for i in score]
                        # # # Check the prediction performance
                        # prediction_2 = [0 if i == -1 else 1 for i in clf.predict(X_test)]
                        # prediction_1 = [0 if i < score_threshold else 1 for i in score]

                        prediction = clf.predict(X_test)
                        print("prediction_2:", Counter(clf.predict(X_test)))
                        prediction[prediction == -1] = 0
                        prediction[prediction == 1] = 1

                        print("true:", Counter(y_test))

                        print('__________________________Results__________________________')
                        # print('tn fp fn tp')
                        # print("__________________________option 1__threshold________________________")
                        # print(confusion_matrix(y_test, prediction_1).ravel())
                        # print("__________________________option 2___________________________")
                        # print(confusion_matrix(y_test, prediction_2).ravel())
                        score = f1_score(y_test, prediction)
                        print('F1 Score: %.3f' % score)
                        precision = precision_score(y_test, prediction, average='binary')
                        print('Precision: %.3f' % precision)
                        recall = recall_score(y_test, prediction, average='binary')
                        print('Recall: %.3f' % recall)
                        # model_confuse_matrix(X_test,y_test, clf, method, clf_dataset, test_dataset,name_classifier)



                    else:
                        prediction = clf.predict(X_test)
                        # shap graph
                        model_shap_plot(X_test, clf, method, clf_dataset, test_dataset,name_classifier, dependence_feature=None)
                        # features importance graph
                        plot_feature_importances(clf, X_test,clf_dataset, test_dataset, name_classifier)
                        model_confuse_matrix(X_test,y_test, clf, method, clf_dataset, test_dataset,name_classifier)

                    if name_classifier == 'isolation_forest':
                        # features shap graph
                        model_shap_plot(X_test, clf, method, clf_dataset, test_dataset,name_classifier, dependence_feature=None)

                   #### ROC Curve ####
                    # y_score = clf.score_samples(X_test)
                    #
                    # from sklearn.metrics import roc_curve
                    # fpr, tpr, thresholds = roc_curve(y_test, y_score)
                    # import matplotlib.pyplot as plt
                    # plt.plot(fpr, tpr, 'k-', lw=2)
                    # plt.xlabel('FPR')
                    # plt.ylabel('TPR')
                    # plt.show()
                    ###################################################################

                    # precision recall curve
                    from sklearn.metrics import precision_recall_curve
                    # from sklearn.metrics import f1_score
                    from sklearn.metrics import auc
                    from matplotlib import pyplot

                    # predict probabilities
                    # lr_probs = clf.decision_function(X_test)
                    # # keep probabilities for the positive outcome only
                    # lr_probs = lr_probs[:, 1]
                    # # predict class values
                    # yhat = clf.predict(X_test)
                    # lr_precision, lr_recall, _ = precision_recall_curve(y_test, lr_probs)
                    # lr_f1, lr_auc = f1_score(y_test, yhat), auc(lr_recall, lr_precision)
                    # # summarize scores
                    # print('Logistic: f1=%.3f auc=%.3f' % (lr_f1, lr_auc))
                    # # plot the precision-recall curves
                    # no_skill = len(y_test[y_test == 1]) / len(y_test)
                    # pyplot.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
                    # pyplot.plot(lr_recall, lr_precision, marker='.', label='Logistic')
                    # # axis labels
                    # pyplot.xlabel('Recall')
                    # pyplot.ylabel('Precision')
                    # # show the legend
                    # pyplot.legend()
                    # # show the plot
                    # pyplot.show()


                    ms = measurement(y_test, prediction)
                    print("TNR:", ms['TNR'])
                    # print(confusion_matrix(y_test, prediction).ravel())
                    # print(Counter(prediction))
                    # print(classification_report(y_test, prediction))

                    if ms_table is None:
                        ms_table = pd.DataFrame(columns=list(ms.keys()), dtype=object)
                    name_method = "model:" + clf_dataset.split("negative")[0] + "/" + "test:" + test_dataset.split("negative")[0]
                    ms_table.loc[name_method] = ms



            except NoModelFound:
                    pass

    res_table.sort_index(axis=0, inplace=True)
    res_table.sort_index(axis=1, inplace=True)

    print(res_table)
    # to_csv(res_table, results_dir / "summary" / "diff_summary_stratify.csv")
    # print(name_classfier)
    to_csv(ms_table, results_dir /"results_iterations" / name_classifier /f"measurement_summary_{number_iteration}.csv")
    print("END result test")
    return ms_table


# different_results_summary(method_split="underSampling", model_dir="models_underSampling")

# different_results_summary(method_split="one_class_svm", model_dir="models_one_class_svm")
