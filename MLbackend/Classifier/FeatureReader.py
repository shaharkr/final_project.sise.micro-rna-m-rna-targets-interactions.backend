from pathlib import Path
from typing import List
from numpy import ndarray
from pandas import DataFrame
from utilsfile import read_csv, to_csv
# from utils.utilsfile import read_csv, to_csv
import pandas as pd

class FeatureReader:
    def __init__(self, expected_num_of_features: int):
        self.expected_num_of_features = expected_num_of_features

    def _feature(self, data) -> List:

        col_list = list(data.columns)
        all_features = col_list[col_list.index("Seed_match_interactions_all"):]
        # all_features = col_list[col_list.index("Acc_P10_10th"):]

        # In this step we read all featuers include hot encoding and after that we remove them
        # in the next step
        assert len(all_features) == 590, f"All feature read error. Read: {len(all_features)}"
        return all_features

    def file_reader(self, in_file: Path) -> (DataFrame, ndarray):
        data: DataFrame = read_csv(in_file)

        # we care on object type columns
        for i in range(1, 21):
            name_col = "miRNAMatchPosition_" + str(i)
            data[name_col] = data[name_col].astype('category').cat.codes
        data["microRNA_name"] = data[name_col].astype('category').cat.codes
        data.drop(columns=['microRNA_name'], inplace=True)

        return self.df_reader(data)

    def df_reader(self, in_df: DataFrame) -> (DataFrame, ndarray):
        y: ndarray = in_df.Label.ravel()
        feature_list: List = self._feature(in_df)
        X = in_df[feature_list]


        assert len(X.columns) == self.expected_num_of_features, f"""Read error. Wrong number of features.
               Read: {len(X.columns)}
               Expected: {self.expected_num_of_features}"""
        return X, y


class AllFeatursReader(FeatureReader):
    def __init__(self):
        super().__init__(580)

    def _feature(self, data) -> List:
        all_features = super()._feature(data)
        return all_features


class HotEncodingReader(FeatureReader):
    def __init__(self):
        super().__init__(90)

    def _feature(self, data) -> List:
        all_features = super()._feature(data)
        return [f for f in all_features if str(f).startswith("HotPairing")]


class AllWithoutEncodingReader(FeatureReader):
    def __init__(self):
        #super().__init__(580-90)
        super().__init__(500)
    def _feature(self, data) -> List:
        all_features = super()._feature(data)
        return [f for f in all_features if not str(f).startswith("HotPairing")]


reader_dict = {
    "all": AllFeatursReader(),
    "hot_encoding": HotEncodingReader(),
    "without_hot_encoding": AllWithoutEncodingReader()
}

reader_selection_parameter = None


def get_reader() -> FeatureReader:
    assert reader_selection_parameter is not None, "reader selection parameter is none"
    return reader_dict[reader_selection_parameter]