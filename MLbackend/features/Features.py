from abc import ABC, abstractmethod
from functools import reduce
from typing import Dict, List, Tuple
from pandas import Series
from duplex.Duplex import Duplex
import pandas as pd
from utils.utilsfile import get_subsequence_by_coordinates


def extract_site_coordinates(site: str, region_sequence: str, start: int, end: int) -> Tuple[int, int]:
    site_len = len(site)

    clean_site = site.replace("*", "").replace("#", "")
    first_valid_nt = [site.find(base) for base in ["A", "U", "C", "G"]]
    first_valid_nt = int(min([i for i in first_valid_nt if i!=-1]))

    target_mrna = get_subsequence_by_coordinates(region_sequence, start, end)
    # print("site:", site)
    # print("target:", target_mrna)
    # print("clean site:", clean_site)
    assert target_mrna.find(clean_site) != -1, f"""cant find the site in full_mrna sequence
    site: {site}
    mrna: {region_sequence}
    start: {start}
    end: {end}
    direct: {region_sequence.find(site)}"""

    new_start = start - first_valid_nt + target_mrna.find(clean_site)
    new_end = new_start + site_len - 1
    return new_start, new_end


class Features(ABC):
    def __init__(self, duplex: Duplex, miRNA_sequence: str, site: str,
                 start: int, end: int, full_mrna: str):
        self._duplex = duplex
        self._miRNA_sequence = miRNA_sequence

        # This not the original site in both the case
        self._site = site

        self._start, self._end = extract_site_coordinates(duplex.site[::-1], full_mrna, int(start), int(end))
        # print(self._start == start)
        # print(self._end == end)

        #print("****************************************************************")
        self._region_sequence = full_mrna

        self._features_dict: Dict = {}

    def get_features(self) -> Series:
        def u(x, y):
            x.update(y)
            return x

        self.extract_features()
        united_dicts = reduce(u, self._features_dict.values())
        return pd.Series(united_dicts)

    @abstractmethod
    def extract_features(self):
        pass

    # site_with_stars_and_hashtags, mRNA_seq_extended, mRNA_seq_extended_offset, full_mrna_seq, ensg,
    # hint)
#
#
# def feature_extraction(miRNA: str, site: str, start: int, end: int, sequence: str, feature_cls_list: List[Features]) -> Series:
#     dp: Duplex = Duplex.fromChimera(miRNA, site)
#     features = [feature_cls(dp, miRNA, site, start, end, sequence).get_features() for feature_cls in feature_cls_list]
#     return pd.concat(features)
#
