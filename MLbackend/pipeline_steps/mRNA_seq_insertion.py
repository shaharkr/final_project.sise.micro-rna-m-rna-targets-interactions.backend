from collections import namedtuple
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame, Series
import click

from consts.mirna_utils import MIRBASE_FILE
from consts.pipeline_steps import READ_PATH, MIRNA_SEQ_PATH, SITE_PATH
from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, read_csv, to_csv


def read_rna_from_fasta(fasta_file: Path) -> DataFrame:
    df = pd.DataFrame(columns=["mRNA ID", "mRNA sequence"])
    with fasta_file.open()  as fasta:
        logger.info(f"read_rna_from_fasta: read fasta file {fasta_file}")
        cnt = 0
        for seq_record in SeqIO.parse(fasta, 'fasta'):  # (generator)
            cnt += 1
            df = df.append(
                pd.Series([seq_record.id, str(seq_record.seq)], index=df.columns), ignore_index=True)
    logger.info(f"read {cnt} rna sequences")
    return df


def rna_insertion(fin_full_path: Path, fout_full_path: Path, rna_df: DataFrame):
    logger.info(f"Insert rna sequence to {fin_full_path}")
    df: DataFrame = read_csv(fin_full_path)
    join_df = df.merge(rna_df, how="left", left_on="mRNA ID", right_on="mRNA ID", validate="many_to_one")
    to_csv(join_df, fout_full_path)
    logger.info(f"Finish the rna sequence insertion to {fin_full_path}")


@click.command()
@click.argument('fasta_file', type=str)
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def insert_from_fasta(fasta_file: str, fin: str, fout:str):
    rna_df: DataFrame = read_rna_from_fasta(Path(fasta_file))
    rna_insertion(Path(fin), Path(fout), rna_df)


@click.group()
def cli():
    pass


cli.add_command(insert_from_fasta)



if __name__ == '__main__':
    cli()

