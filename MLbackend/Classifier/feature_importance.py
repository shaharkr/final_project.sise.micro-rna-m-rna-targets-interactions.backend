import pickle
from pathlib import Path
from consts.global_consts import ROOT_PATH, NEGATIVE_DATA_PATH, MERGE_DATA, DATA_PATH_INTERACTIONS
from utilsfile import read_csv, to_csv
import xlsxwriter
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import FunctionTransformer
import Classifier.FeatureReader as FeatureReader
from Classifier.FeatureReader import get_reader
from utilsClassifier import mean_std


DATASET_LIST = ['human_dataset1', 'human_dataset2', 'human_dataset3', 'mouse_dataset1',  'mouse_dataset2',
                'celegans_dataset1', 'celegans_dataset2', 'cattle_dataset1']
dataset_translation_list = {'human_dataset1': "h1", 'human_dataset2': "h2", 'human_dataset3': "h3",
                                    'mouse_dataset1': "m1", 'mouse_dataset2': "m2",
                                    'celegans_dataset1': "ce1", 'celegans_dataset2': "ce2",
                                    'cattle_dataset1': "ca1"}

DATASET_LIST.sort()

def generate_importance_files(method_split, model_dir):
    results_dir =  ROOT_PATH / Path("Results")
    method = 'xgbs'
    results_dir_models = ROOT_PATH / Path("Results") / model_dir

    importance_types = ["weight", "gain", "cover", "total_gain", "total_cover"]
    test_dir = DATA_PATH_INTERACTIONS / "test" / method_split
    FeatureReader.reader_selection_parameter = "without_hot_encoding"
    feature_reader = get_reader()
    for f_test in test_dir.glob("*test*"):
        f_stem = f_test.stem
        test_dataset = f_stem.split(".csv")[0]
        X_test, y_test = feature_reader.file_reader(test_dir / f"{test_dataset}.csv")
        features= list(X_test.columns)

    clf_datasets = [f.stem for f in results_dir_models.glob("*_xgbs*model")]

    for model_file in clf_datasets:
        print(f"model file = {model_file}")
        feature_importance = pd.DataFrame(index=features)
        model_file= model_file+".model"
        clf_file = results_dir / model_dir/model_file

        with clf_file.open("rb") as f:
            clf = pickle.load(f)
        for it in importance_types:
            feature_importance[it] = pd.Series(clf.get_booster().get_score(importance_type=it))
        feature_importance.fillna(value=0, inplace=True)
        to_csv(feature_importance, results_dir / "feature_importance" / f"feature_importance_{model_file}.csv")
        print("save feature importance file")


def feature_importance_table():
    writer = pd.ExcelWriter('Results/feature_importance.xlsx', engine='xlsxwriter') # supplementary_file
    for dataset in DATASET_LIST:
        print(dataset)
        results_dir = ROOT_PATH / Path("Results")

        mstd_df, m, _ = mean_std(dir=results_dir, match="self*", summary_shape=(502, 5),
                                 summary_file_name=f"feature_importance_{dataset}_xgbs_no_encoding.csv")
        m.to_csv(Path("Results")/ f"new_feature_importance_{dataset}.csv")
        mstd_df.sort_index().to_excel(writer, sheet_name=dataset_translation_list[dataset])
    writer.save()


def create_feature_df():
    dir = Path("Results")
    fi = pd.DataFrame()
    for f in dir.glob("*new*"):
        dataset = dataset_translation_list[f.stem.split("importance_")[1]]
        fi[dataset] = pd.read_csv(f, index_col=0)["gain"]
    return fi


def feature_importance_plot(fi):
    fi_sort = pd.DataFrame(columns=fi.columns)
    for c in fi_sort.columns:
        fi_sort[c] = fi[c].sort_values(ascending=False).values

    # Plot the full df
    fi_sort.plot()
    plt.xlabel('Features')
    plt.ylabel('Gain')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.savefig(dir/ "feature_importance_full.pdf", format="pdf", bbox_inches='tight')

    # Plot the zoom df
    zoom_df = fi_sort.head(20)
    zoom_df.plot()
    plt.xlabel('Features')
    plt.ylabel('Gain')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.xticks(range(0, len(zoom_df), 5))
    plt.savefig(dir / "feature_importance_zoom.pdf", format="pdf", bbox_inches='tight')


def max_transformer(X):
    scale = 100 / X.max(axis=0)
    return scale * X

def get_top_features (fi, n):
    top_features = set()
    for c in fi.columns:
        current_dataset_top_features = set(fi[c].sort_values(ascending=False).head(n).index)
        top_features = top_features.union(current_dataset_top_features)
    top_df = fi.loc[top_features,:]
    transformer = FunctionTransformer(max_transformer)
    top_df = transformer.transform(top_df)
    top_df["mean"] = top_df.mean(axis=1)
    top_df.sort_values(by="mean", ascending=False, inplace=True)
    top_df = top_df.round(0).astype(int)

    top_df.to_csv("Results/top_features.csv")


def main():
    importance_types = ["weight", "gain", "cover", "total_gain", "total_cover"]
    # generate_importance_files(method_split="stratify", model_dir="models_stratify")
    feature_importance_table()
    fi = create_feature_df()
    feature_importance_plot(fi)
    get_top_features(fi, 6)




if __name__ == '__main__':
    main()
