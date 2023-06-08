from abc import ABC
from itertools import product
from typing import Dict

from pandas import Series

from consts.features import HOT_ENCODING_LEN
from features.Features import Features
from utils.utilsfile import get_subsequence_by_coordinates


def seq_composition(s: str, prefix: str) -> Dict[str, float]:

    # make change
    # if len(s) < 1:
    #     return {}
    # monomer
    keys = ["".join(p) for p in product('ACGU', repeat=1)]
    monomer_dict = dict.fromkeys(keys, 0)
    for i in range(len(s)):
        k = s[i]
        monomer_dict[k]+=1
    monomer = sum(monomer_dict.values())

    # if len(s) < 2:
    #     return {}
    #dimer
    keys = ["".join(p) for p in product('ACGU', repeat=2)]
    if s.find("#")== True:
        print("f")
    dimer_dict = dict.fromkeys(keys, 0)
    for i in range(len(s) - 1):
        k = s[i:i+2]
        dimer_dict[k]+=1
    dimer = sum(dimer_dict.values())

    # make change
    if monomer == 0:
        monomer = 1
    if dimer == 0:
        dimer = 1
    result = {f"{prefix}_{k}_comp": round(float(v) / monomer, 4) for k, v in monomer_dict.items()}
    result.update({f"{prefix}_{k}_comp": round(float(v) / dimer, 4) for k, v in dimer_dict.items()})
    return result


class MrnaFeatures(Features):
    def extract_features(self):
        self._features_dict["dte"] = self.distance_to_end()
        self._features_dict["dts"] = self.distance_to_start()
        # make change -
        # Gilad
        self._features_dict["target_composition"] = seq_composition(self._site, "MRNA_Target")
        # Efrat
        # self._features_dict["target_composition"] = seq_composition(self._duplex.site[::-1], "MRNA_Target")

        try:
            # make change
            mrna_up = get_subsequence_by_coordinates(self._region_sequence, max(self._start - 71, 0), self._start - 1)
            self._features_dict["flanking_up_composition"] = seq_composition(mrna_up, "MRNA_Up")
        except ValueError:
            # self._features_dict["flanking_up_composition"] = {}
            self._features_dict["flanking_up_composition"] = seq_composition("", "MRNA_Up")

        try:
            # make change
            mrna_down = get_subsequence_by_coordinates(self._region_sequence, self._end + 1, min(self._end + 71,len(self._region_sequence)))
            if len(mrna_down) == 1:
                print("ffff")
            self._features_dict["flanking_down_composition"] = seq_composition(mrna_down, "MRNA_Down")
        except ValueError:
            # self._features_dict["flanking_down_composition"] = {}
            self._features_dict["flanking_down_composition"] = seq_composition("", "MRNA_Down")




        self._features_dict["hot_encoding"] = self.pair_hot_encoding(HOT_ENCODING_LEN)

    def distance_to_end(self) -> Dict[str, float]:
        mrna_len = len(self._region_sequence)
        return {'MRNA_Dist_to_end': round(float(mrna_len - self._end - 1) / mrna_len, 4)}

    def distance_to_start(self) -> Dict[str, float]:
        mrna_len = len(self._region_sequence)
        return {'MRNA_Dist_to_start': round(float(self._start) / mrna_len, 4)}


    # # 8. pair hot-encoding
    def pair_hot_encoding(self, length: int):
        mirna = self._miRNA_sequence[:length]
        # if we need- bug! duplex.site
        mrna = self._site[:length]

        def hot_coding(seq):
            if seq == 'A' or seq == 'a':
                he = [1, 0, 0, 0, 0]
            elif seq == 'U' or seq == 'u':
                he = [0, 1, 0, 0, 0]
            elif seq == 'T' or seq == 't':
                he = [0, 1, 0, 0, 0]
            elif seq == 'G' or seq == 'g':
                he = [0, 0, 1, 0, 0]
            elif seq == 'C' or seq == 'c':
                he = [0, 0, 0, 1, 0]
            else:
                he = [0, 0, 0, 0, 1]
            return he

        PHE = {}
        for i in range(len(mirna)):
            for j in range(5):
                key = 'HotPairingMirna_he_P%s_L%s' % (str(i + 1), str(j + 1))
                PHE[key] = hot_coding(mirna[i])[j]

        for i in range(len(mrna)):
            for j in range(5):
                key = 'HotPairingMRNA_he_P%s_L%s' % (str(i + 1), str(j + 1))
                PHE[key] = hot_coding(mrna[i])[j]
        return PHE
