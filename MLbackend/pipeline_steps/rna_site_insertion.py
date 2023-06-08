from collections import namedtuple
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame, Series
import click

from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, SITE_EXTRA_CHARS
from consts.mirna_utils import MIRBASE_FILE
from consts.pipeline_steps import READ_PATH, MIRNA_SEQ_PATH, SITE_PATH
from utils.genome import extract_seq_from_chromosome
from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, get_wrapper, read_csv, to_csv


@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def insert_site_by_coordinates(fin: str, fout:str):
    logger.info(f"Insert site to {fin}")
    df: DataFrame = read_csv(Path(fin))
    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "mRNA sequence", "chimera_start", "chimera_end",
                                           extra_chars=SITE_EXTRA_CHARS),
                          axis=1)

    to_csv(df, Path(fout))
    logger.info(f"finish the site sequence insertion to {fin}")

@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def get_site_from_extended_site(fin: str, fout:str):
    def calc_chimera_start(seq: str, subseq: str) -> int:
        try:
            if seq.find(subseq) == -1:
                return -1
            return seq.find(subseq) + 1
        except AttributeError:
            return -1

    def calc_chimera_end(chimera_start: int, seq_extended: str) -> int:
        if chimera_start == -1:
            return -1
        return chimera_start + len(seq_extended) - 1 - HUMAN_SITE_EXTENDED_LEN

    logger.info(f"Insert site to {fin}")
    df: DataFrame = read_csv(Path(fin))
    df["chimera_start"] = df.apply(func=get_wrapper(calc_chimera_start,
                                                    'region sequence', 'mRNA_seq_extended'),
                                   axis=1)
    df["chimera_end"] = df.apply(func=get_wrapper(calc_chimera_end,
                                                  'chimera_start', 'mRNA_seq_extended'),
                                 axis=1)

    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "region sequence", "chimera_start", "chimera_end",
                                           extra_chars=SITE_EXTRA_CHARS),
                          axis=1)

    to_csv(df, Path(fout))
    logger.info(f"finish the site sequence insertion to {fin}")


@click.argument('fin', type=str)
def insert_site_from_chromosome(fin: str, fout: str, chr_dir: str):
    logger.info(f"Insert site from chromosome to {fin}")
    df: DataFrame = read_csv(Path(fin))
    df["site"] = df.apply(func=get_wrapper(extract_seq_from_chromosome,
                                           'chr', 'start', 'end', 'strand',
                                           directory=Path(chr_dir)),
                          axis=1)
    df["site"] = df["site"].apply(lambda x: x.upper())

    to_csv(df, Path(fout))
    logger.info(f"finish the site sequence insertion to {fin}")



@click.group()
def cli():
    pass


cli.add_command(insert_site_by_coordinates)
cli.add_command(get_site_from_extended_site)
cli.add_command(insert_site_from_chromosome)


if __name__ == '__main__':
    cli()

