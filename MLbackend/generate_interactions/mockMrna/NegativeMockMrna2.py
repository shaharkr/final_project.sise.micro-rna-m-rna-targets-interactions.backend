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
from ushuffle import shuffle, Shuffler


""" In this method, we generate mockMrna(by shuffle the +-100 around the original site.
  Then, we calculate the duplex of the mockMrna with the original mirna"""


class MockMRNA(object):

    def __init__(self, organism, tmp_dir, min_num_of_pairs):
        self.tmp_dir = tmp_dir
        self.min_num_of_pairs = min_num_of_pairs

    def generate_mrna_mock_denucleotides(self, mrna, th=5):

        seq = list(mrna.replace('T', 'U').upper())
        seq = ''.join(seq)
        seq = seq.encode('utf-8')

        num_shuffle = 0
        equal_to_itself = True
        while equal_to_itself:
            # seq_couple = [seq[x:x + 2] for x in range(0, len(seq), 2)]
            # random.shuffle(seq_couple)
            shuffler = Shuffler(seq, 2)
            seq_byte = shuffler.shuffle()
            seq_str = seq_byte.decode("utf-8")
            num_shuffle += 1
            if num_shuffle % 10000 == 0:
                print(num_shuffle)
            if num_shuffle > 100000:
                break
            # flat_seq_mock = [item for sublist in seq_couple for item in sublist]
            equal_to_itself = seq_str == seq

        mrna_mock = ''.join(seq_str)
        return mrna_mock

    def generate_mrna_mock_nucleotides_ushuffle(self, mrna, th=5):

        seq = list(mrna.replace('T', 'U').upper())
        seq = ''.join(seq)
        seq = seq.encode('utf-8')

        num_shuffle = 0
        equal_to_itself = True
        while equal_to_itself:
            # seq_couple = [seq[x:x + 2] for x in range(0, len(seq), 2)]
            # random.shuffle(seq_couple)
            shuffler = Shuffler(seq, 1)
            seq_byte = shuffler.shuffle()
            seq_str = seq_byte.decode("utf-8")
            num_shuffle += 1
            if num_shuffle % 10000 == 0:
                print(num_shuffle)
            if num_shuffle > 100000:
                break
            # flat_seq_mock = [item for sublist in seq_couple for item in sublist]
            equal_to_itself = seq_str == seq

        mrna_mock = ''.join(seq_str)
        return mrna_mock

    # def generate_mrna_mock_nucleotides(self, mrna, th=5):
    #
    #     seq = list(mrna.replace('T', 'U').upper())
    #     seq_original = list(mrna.replace('T', 'U').upper())
    #
    #     num_shuffle = 0
    #     equal_to_itself = True
    #     while equal_to_itself:
    #         random.shuffle(seq)
    #         num_shuffle += 1
    #         if num_shuffle % 10000 == 0:
    #             print(num_shuffle)
    #         if num_shuffle > 100000:
    #             break
    #         equal_to_itself = seq_original == seq
    #
    #     mrna_mock = ''.join(seq)
    #     return mrna_mock
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

    def insert_mock_site(self, full_mrna, start, end, site):

        full_mrna = full_mrna[:start] + site + full_mrna[end + 1:]
        return full_mrna

    def generate_negative_seq(self, original_mirna, full_mrna_original, site, start, end, name_shuffle, num_of_tries=10000):

        for i in range(num_of_tries):

            # Take the +- 100 around the site
            site_surrounding100 = get_subsequence_by_coordinates(full_mrna_original, start, end, extra_chars=100)
            new_start, new_end = get_substring_index(full_mrna_original, site_surrounding100)

            # Shuffle the sequence
            if name_shuffle == 'nucleotides':
                mock_sub_mrna = self.generate_mrna_mock_nucleotides_ushuffle(site_surrounding100)
            else:
                mock_sub_mrna = self.generate_mrna_mock_denucleotides(site_surrounding100)

            # Insert the mock site to the original mrna
            full_mrna = self.insert_mock_site(full_mrna_original, new_start, new_end, mock_sub_mrna)

            # find new site on mrna_surrounding50
            mrna_site: str = get_subsequence_by_coordinates(full_mrna, start, end, extra_chars=50)

            # this step check that the mock is vaild. The meaning is that the interaction is canonical
            # or noncanonical

            canonic_seed, non_canonic_seed, num_of_pairs = self.valid_negative_seq(original_mirna, mrna_site)
            cond1 = canonic_seed
            cond2 = non_canonic_seed
            if cond1 or cond2:
                properties = {
                    "mirna": original_mirna,
                    "full_mrna": full_mrna,
                    "canonic_seed": canonic_seed,
                    "non_canonic_seed": non_canonic_seed,
                    "num_of_pairs": num_of_pairs,
                    'site':  mrna_site
                }
                return True, properties
        return False, {}


def worker(organism, fin, fout_name, tmp_dir,name_shuffle):

    print("##################NEW FILE#################################")
    print(fin)
    fout = filename_suffix_append(fout_name, "_negative")
    min_num_of_pairs = CONFIG["minimum_pairs_for_interaction"]
    ns = MockMRNA(organism, tmp_dir=tmp_dir, min_num_of_pairs=min_num_of_pairs)
    in_df = read_csv(fin)
    in_df = in_df[(in_df["Seed_match_canonical"] == 'True') | (in_df["Seed_match_noncanonical"] == 'True')]
    neg_df = pd.DataFrame()

    i=0
    in_df.rename(columns={'miRNA ID': 'miRNA_ID'}, inplace=True)

    for index, row in in_df.iterrows():
        print(f"$$$$$$$$$$$$$$$ {i} $$$$$$$$$$$$$$$$$$4")
        i += 1
        valid, properties = ns.generate_negative_seq(row['miRNA sequence'], row['sequence'],
                                                     row['site'], row['start'], row['end'], name_shuffle)

        if not valid:
            continue

        new_row = pd.Series()
        new_row['paper name'] = row['paper name']
        new_row['organism'] = row['organism']
        new_row['miRNA ID'] = row.miRNA_ID
        new_row['miRNA sequence'] = row['miRNA sequence']
        new_row['Gene_ID'] = row.Gene_ID
        new_row['full_mrna'] = properties['full_mrna']
        new_row["num_of_pairs"] = properties["num_of_pairs"]
        new_row["site"] = properties["site"]

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


def main(name_shuffle):
    file_name = ROOT_PATH / "data/positive_interactions/positive_interactions_merge"
    tmp_base = ROOT_PATH / "generate_interactions/mockMrna/"
    print("tmp:", tmp_base)
    files = list(file_name.glob('**/*.csv'))
    for p in files:
        fout_name = p.name.split('.csv')[0] + "_" + name_shuffle + '_method2.csv'
        if "darnell_human" in fout_name:
            worker('hsa', p, fout_name, tmp_base, name_shuffle)


# main()




