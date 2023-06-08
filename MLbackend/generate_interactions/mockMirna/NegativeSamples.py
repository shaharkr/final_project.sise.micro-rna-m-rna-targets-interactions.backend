from pathlib import Path
from datetime import datetime
from utils.logger import logger
from consts.mirna_utils import MIRBASE_FILE
from consts.global_consts import DUPLEX_DICT
from utils.utilsfile import read_csv, to_csv
import pandas as pd
import numpy as np
from utils.utilsfile import *
from consts.global_consts import ROOT_PATH, DATA_PATH
from duplex.ViennaDuplex import ViennaDuplex
from duplex.ViennaDuplex import *
import random
from features.SeedFeatures import *
#import MirBaseUtils.mirBaseUtils as MBU
import mirna_utils.mirbase as MBU
from multiprocessing import Process
from consts.global_consts import CONFIG
from duplex.Duplex import Duplex


class NegativeSamples(object):

    def __init__(self, organism, tmp_dir, min_num_of_pairs):
        self.tmp_dir = tmp_dir
        self.min_num_of_pairs = min_num_of_pairs
        mirna_prefix_table = {"human" : "hsa",
                               "mouse" : "mmu",
                               "elegans" : "cel",
                               "celegans": "cel",
                               "cattle": "bta"}

        #miRBase_dic = MBU.read_mirbase_file(mirna_prefix_table[organism.lower()])
        file_name = MIRBASE_FILE
        logger.info(f"Reading file {file_name}")
        df_mirna = pd.read_csv(file_name)
        df_mirna_hsa = df_mirna[df_mirna['miRNA ID'].str.contains('hsa')]
        self.miRBase_seq_list = df_mirna_hsa["miRNA sequence"].tolist()

    # # 4. a function used to generate mock mirna which their seed never equal to any mirna in the miRBase.

    def generate_mirna_mock(self, mirna, th=5):
        def seed_equal(a, b):
            return sum(a[i] == b[i] for i in range(len(a)))

        def equal_to_mirbase(mir, th=5):
            for seq in self.miRBase_seq_list:
                seq = list(seq)
                e27 = seed_equal(mir[1:7], seq[1:7])
                if e27 > th:
                    return True
                e38 = seed_equal(mir[2:8], seq[2:8])
                if e38 > th:
                    return True
            return False

        mirna_list_o = list(mirna.replace('T', 'U').upper())
        mirna_list_r = list(mirna.replace('T', 'U').upper())
        equal_to_itself = True

        num_shuffle = 0
        while equal_to_itself or eq_mirbase:
            random.shuffle(mirna_list_r)
            num_shuffle += 1
            if num_shuffle % 10000 == 0:
                print(num_shuffle)
            if num_shuffle > 100000:
                break

            # check if it equals to itself
            e27 = seed_equal(mirna_list_r[1:7], mirna_list_o[1:7])
            e38 = seed_equal(mirna_list_r[2:8], mirna_list_o[2:8])
            equal_to_itself = e27 > th or e38 > th
            # check against mirbase
            eq_mirbase = equal_to_mirbase(mirna_list_r)

        mirna_m = ''.join(mirna_list_r)
        return mirna_m

    def valid_negative_seq(self, mir, mrna):

        duplex_cls: Duplex = DUPLEX_DICT['ViennaDuplex']
        logger.info(f"{ViennaDuplex} do_duplex")
        dp = duplex_cls.fromChimera(mir, mrna)
        try:
            canonic_seed = dp.canonical_seed
            non_canonic_seed = dp.noncanonical_seed

        except SeedException:
            canonic_seed = False
            non_canonic_seed = False

        # warning: number of pair was before to interactions count
        return canonic_seed, non_canonic_seed, dp.interaction_count

    def generate_negative_seq(self, orig_mirna, full_mrna, num_of_tries=10000):

        for i in range(num_of_tries):
            # check that the mock mirna doesn't appear in mirbase list ot of exists mirna
            mock_mirna = self.generate_mirna_mock(orig_mirna)

            # this step check that the mock is valid. The meaning is that the interaction is canonical
            # or noncanonical

            canonic_seed, non_canonic_seed, num_of_pairs = self.valid_negative_seq(mock_mirna, full_mrna)
            cond1 = canonic_seed
            cond2 = non_canonic_seed
            if cond1 or cond2:
                properties = {
                    "mock_mirna": mock_mirna,
                    "full_mrna": full_mrna,
                    "canonic_seed": canonic_seed,
                    "non_canonic_seed": non_canonic_seed,
                    "num_of_pairs": num_of_pairs
                }
                return True, properties
        return



def worker(organism, fin, fout_name, tmp_dir):

    print("##################NEW FILE#################################")
    print(fin)
    fout = filename_suffix_append(fout_name, "_negative")
    min_num_of_pairs = CONFIG["minimum_pairs_for_interaction"]
    ns = NegativeSamples(organism, tmp_dir=tmp_dir, min_num_of_pairs=min_num_of_pairs)
    in_df = read_csv(fin)
    in_df = in_df[(in_df["Seed_match_canonical"] == 'True') | (in_df["Seed_match_noncanonical"] == 'True')]
    neg_df = pd.DataFrame()

    i=0
    in_df.rename(columns={'miRNA ID': 'miRNA_ID'}, inplace=True)

    for index, row in in_df.iterrows():
        print(f"$$$$$$$$$$$$$$$ {i} $$$$$$$$$$$$$$$$$$4")
        i += 1
        valid, properties = ns.generate_negative_seq(row['miRNA sequence'], row['sequence'])
        if not valid:
            continue

        new_row = pd.Series()
        new_row['paper name'] = row['paper name']
        new_row['organism'] = row['organism']
        new_row['miRNA ID'] = "mock " + row.miRNA_ID
        new_row['miRNA sequence'] = properties["mock_mirna"]
        #new_row['target sequence'] = row['sequence']
        #new_row['number of reads'] = row['number of reads']
        # new_row["canonic_seed"] = properties["canonic_seed"]
        # new_row["non_canonic_seed"] = properties["non_canonic_seed"]
        new_row['Gene_ID'] = row.Gene_ID
        # new_row['mRNA_start'] = 0
        # new_row['mRNA_end'] = len(row['sequence'])
        new_row['full_mrna'] = row['sequence']
        new_row["num_of_pairs"] = properties["num_of_pairs"]

        neg_df = neg_df.append(new_row, ignore_index=True)
        # if i > 10:
        #     break


    ########################
    # Save df to CSV
    ########################
    neg_df.reset_index(drop=True, inplace=True)
    #neg_df.drop(labels=['Seed_match_noncanonical', 'Seed_match_canonical'], axis=1, inplace=True)
    drop_unnamed_col(neg_df)
    neg_df["key"] = neg_df.reset_index().index
    fout = tmp_dir / fout
    to_csv(neg_df, fout)
    print("save:", fout)


def main():
    #INPUT DIR
    file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_merge"
    #OUTPUT
    tmp_base = ROOT_PATH / "generate_interactions/mockMirna/"
    print("tmp:", tmp_base)
    files = list(file_name.glob('**/*.csv'))
    for p in files:
        name = Path("/sise/home/efrco/efrco-master/data/positive_interactions/positive_interactions_merge/darnell_human_ViennaDuplex_features.csv")
        if p != name:
            continue
        fout_name = p.name.split('.csv')[0] + '.csv'
        worker('hsa', p, fout_name, tmp_base)

# main()










