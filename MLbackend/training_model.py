import os
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix
import pickle

def prepare_training_data(path_to_positive, path_to_negative, ratio=0.8,features_to_remove=None):
    """
    :param path_to_positive: Path to positive samples of the miRNA sequences after feature extraction
    :param path_to_negative: Path to negative samples of the miRNA sequences after feature extraction
    :param ratio: ratio of train from the whole data
    :return: two DataFrames: train, test
    """

    # Load the CSV files with the specified data types
    positive_df = pd.read_csv(path_to_positive, low_memory=False)
    negative_df = pd.read_csv(path_to_negative, skiprows=[1])
    if features_to_remove is None:
            features_to_remove = ['miRNA ID', 'miRNA sequence', 'Gene_ID', 'new_key', 'MiRBase ID',
                                  'Node of origin (locus)', 'Node of origin (family)', 'Family Age_Group',
                                  'Family Age', 'Locus Age', 'seed_family', 'site', 'region', 'paper region',
                                  'sequence', 'Gene_name', 'mrna_bulge', 'mrna_inter', 'mir_inter',
                                  'mir_bulge', 'duplex_method'
                                  ]

    # removing features
    features_to_remove_positive = list(set(positive_df.columns).intersection(features_to_remove))
    features_to_remove_negative = list(set(negative_df.columns).intersection(features_to_remove))

    positive_df = positive_df.drop(columns=features_to_remove_positive)
    negative_df = negative_df.drop(columns=features_to_remove_negative)

    # Find common columns in both DataFrames
    common_columns = list(set(positive_df.columns).intersection(negative_df.columns))

    # Remove any columns that are not in both DataFrames
    positive_df = positive_df[common_columns]
    negative_df = negative_df[common_columns]

    # Remove datatypes that are not numeric
    positive_df = positive_df.select_dtypes(include=['float64', 'int64'])
    negative_df = negative_df.select_dtypes(include=['float64', 'int64'])

    # Ensure equal number of positive and negative samples
    num_rows = min(positive_df.shape[0], negative_df.shape[0])
    positive_df = positive_df.head(num_rows)
    negative_df = negative_df.head(num_rows)

    # Add a new column named "target" with 1 for positive and 0 for negative
    positive_df["target"] = 1
    negative_df["target"] = 0

    # Concatenate the DataFrames
    concatenated_df = pd.concat([positive_df, negative_df], ignore_index=True)

    train = concatenated_df.sample(frac=ratio, random_state=42)
    test = concatenated_df.drop(train.index)

    return train, test


def train_model(train_df,
                model_name=None,
                save_path=None,
                save_model=False,
                ):


    X_train = train_df.drop(['target', 'index', 'key'], axis=1)
    y_train = train_df['target']

    # Define the XGBoost model and fit it to the training data
    xgb_model = xgb.XGBClassifier(objective='binary:logistic')
    xgb_model.fit(X_train, y_train)

    if save_model:
        # save the model
        if save_path is None:
            save_path = os.getcwd()
        if model_name is None:
            model_name = 'XGB_model'
        model_file = os.path.join(save_path, f"{model_name}.model")
        with open(model_file, 'wb') as f:
            pickle.dump(xgb_model, f)
        print(f"Model saved to: {model_file}")
    # return the trained model
    return xgb_model


def evaluate_model(model=None, model_path=None, eval_set=None, plot_CM=False):
    # if model is None and model_path is None:
    #     raise ValueError("Either a model or a model_path must be provided.")

    X_eval = eval_set.drop(['target', 'index', 'key'], axis=1)
    y_eval = eval_set['target']
    y_pred = model.predict(X_eval)

    # Calculate the accuracy score
    acc_score = accuracy_score(y_eval, y_pred)
    print(f"Accuracy Score: {acc_score}")

    # Plot the confusion matrix if plot_CM is True
    if plot_CM:
        cm = confusion_matrix(y_eval, y_pred)
        sns.heatmap(cm, annot=True, cmap="Blues", fmt="d")
        plt.title("Confusion Matrix")
        plt.xlabel("Predicted labels")
        plt.ylabel("True labels")
        plt.show()


if __name__ == '__main__':

    path_to_positive = r"C:\Users\user\Desktop\positive_features\pairing_beyond_features.csv"
    path_to_negative = r"C:\Users\user\Desktop\negative_features\mockMirna_pairing_beyond_negative_features.csv"
    train, test = prepare_training_data(path_to_positive, path_to_negative, ratio=0.8)
    print('finished preparing data')
    xgb_model = train_model(train, model_name='Worm',save_model=True)
    print('trained model')
    evaluate_model(model=xgb_model, eval_set=test, plot_CM=True)
