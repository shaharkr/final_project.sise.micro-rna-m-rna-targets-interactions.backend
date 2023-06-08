from abc import ABC
from itertools import product
from typing import Dict

from pandas import Series

from consts.features import HOT_ENCODING_LEN
from duplex.Duplex import Duplex
from duplex.utils import mix_inter_bulge_seq
from features.Features import Features

AU = ['AU', 'UA']
GC = ['GC', 'CG']
GU = ['GU', 'UG']
MM = ['AA', 'AG', 'AC', 'UU', 'UC', 'GA', 'GG', 'CA', 'CU', 'CC']
c = ['AU', 'UA', 'GC', 'CG']
w = ['GU', 'UG']

class SeedException (Exception):
    pass


def startingA(seed: Duplex):
    mirna0 = mix_inter_bulge_seq(seed.mir_inter[0], seed.mir_bulge[0])
    mrna0 = mix_inter_bulge_seq(seed.mrna_inter[0], seed.mrna_bulge[0])

    A_mirna = mirna0 == 'A' and mirna0 + mrna0 not in c
    A_mrna = mrna0 == 'A' and mirna0 + mrna0 not in c
    return 1 if (A_mirna or A_mrna) else 0

def mismatch(seed: Duplex, i):
    GU = ['GU', 'UG']
    pair: str = seed.mir_bulge[i] + seed.mrna_bulge[i]
    # if pair in GU:
    #     return 0
    return 1 if len(pair.strip()) == 2 else 0

def bulge(a, b):
    return sum([a[i] != ' ' and b[i] == ' ' for i in range(len(a))])

def countGU(seed: Duplex):
    return sum([pair in w for pair in seed.pair_iterator()])

def startingIndex(seed: Duplex):
    for i, (j, p) in enumerate(seed.mir_iterator()):
        if seed.mir_inter[j] != ' ':
            return i + 1
    return -1
    raise SeedException("not valid seed. No interaction at all")





class SeedFeatures(Features):
    def __init__(self, duplex: Duplex, miRNA_sequence: str, site: str, start: int, end: int, region_sequence: str):
        super().__init__(duplex, miRNA_sequence, site, start, end, region_sequence)
        self._seed = duplex.seed
        self._seed_2_7 = duplex.extract_seed(2, 7)
        self._seed_3_8 = duplex.extract_seed(3, 8)

    def extract_features(self):
        seed_features_dict = {}

        seed_features_dict['Seed_match_canonical'] = self._duplex.canonical_seed
        seed_features_dict['Seed_match_noncanonical'] = False if self._duplex.canonical_seed else self._duplex.noncanonical_seed

        seed_regions = {"all" : self._seed,
                        "2_7" : self._seed_2_7,
                        "3_8" : self._seed_3_8}
        

        for key, seed in seed_regions.items():
            seed_features_dict[f"Seed_match_interactions_{key}"] = seed.interaction_count
            seed_features_dict[f"Seed_match_GU_{key}"] = countGU(seed)

        seed: Duplex = self._seed

        seed_features_dict['Seed_match_A'] = startingA(seed)
        seed_features_dict['Seed_match_start'] = startingIndex(seed)

        mirna_last_nt = list(seed.mir_iterator())[-1][0]
        seed_features_dict['Seed_match_mismatch_left'] = mismatch(seed, 0)
        seed_features_dict['Seed_match_mismatch_right'] = mismatch(seed, mirna_last_nt)

        s_i = seed.mir_bulge.find(' ')
        s_e = seed.mir_bulge.rfind(' ')
        seed_features_dict['Seed_match_mismatch_inner'] = 0

        for i in range(s_i, s_e + 1):
            seed_features_dict['Seed_match_mismatch_inner'] += mismatch(seed, i)

        seed_features_dict['Seed_match_bulge_target'] = bulge(seed.mrna_bulge, seed.mir_bulge)
        seed_features_dict['Seed_match_bulge_mirna'] = bulge(seed.mir_bulge, seed.mrna_bulge)

        self._features_dict["seed_features"] = seed_features_dict



# 
# #
# # def count_not_space(s):
# #     return len(s) - s.count(' ')
# #
# #
# #
# #         # self.mirna_last_nt = list(seed.mir_iterator())[-1][0]
# #
# #     def extract_seed_features(self):
# #         self.seed_match_type()
# 
# 
#     def valid_seed(self):
#         # Valid Seed must have at least 6 combination of interactions and GUs
#         assert self.smt_dic is not None, "The seed dict hasn't initiated yet."
#         # return (self.smt_dic['Seed_match_compact_interactions'] + self.smt_dic['Seed_match_compact_GU']) >= 6
#         return self.canonical_seed() or self.non_canonical_seed()
