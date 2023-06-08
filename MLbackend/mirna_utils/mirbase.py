from ftplib import FTP
from pathlib import Path
from typing import List
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO
import click

import gzip

from consts.mirna_utils import MIRBASE_URL, VERSIONS_DIR, FASTA_DIR, MIRBASE_FILE, MATURE_FILE
from utils.logger import logger


def mirbase_connect() -> FTP:
    logger.info("Connecting to ftp")
    ftp = FTP(MIRBASE_URL)
    ftp.login()
    ftp.cwd(VERSIONS_DIR)
    return ftp


def get_mirbase_versions(ftp: FTP, min_version: float) -> List[str]:
    def my_float(x):
        try:
            return (float(x))
        except ValueError:
            return None

    df: DataFrame = DataFrame(ftp.nlst(), columns=["key"])
    df["float"] = df["key"].apply(my_float)
    df = df[df["float"] >= min_version]
    df.sort_values(by="float", ascending=False, inplace=True)
    return list(df["key"])


def get_mature_file_from_dir_list(ftp: FTP, dir_list: List[str]):
    logger.info("Browsing the FTP")
    for d in dir_list:
        ftp.cwd(d)
        logger.info(f"Changed to {ftp.pwd()}")
        target_file = FASTA_DIR / f"{MATURE_FILE}_ver{d}"
        ftp.retrbinary("RETR " + MATURE_FILE, open(target_file, 'wb').write)
        ftp.cwd('..')

@click.command()
@click.argument('min_version', type=float)
def read_fasta_files(min_version: float):
    mirbase_df = pd.DataFrame(columns=["version", "miRNA ID", "miRNA sequence"])

    for fasta_gz in FASTA_DIR.glob("*"):

        ver = str(fasta_gz).split("_ver")[1]
        if float(ver) < min_version:
            continue
        logger.info(f"read fasta file {fasta_gz}")

        with gzip.open(fasta_gz, 'rt') as fasta:
            for seq_record in SeqIO.parse(fasta, 'fasta'):  # (generator)
                mirbase_df = mirbase_df.append(
                    pd.Series([ver, seq_record.id, str(seq_record.seq)], index=mirbase_df.columns), ignore_index=True)

    mirbase_df.to_csv(MIRBASE_FILE)
    logger.info("save mirbase file")



@click.command()
def postprocess():
    def get_prefix(x):
        return x.split("-")[0]
    df = pd.read_csv(MIRBASE_FILE, index_col=0)
    #df.sort_values(by="version", ascending=False, inplace=True)
    logger.info(f"mirbase file before postprocess {df.shape}")
    df.drop_duplicates(subset=["miRNA ID"], keep='first', inplace=True)
    df["prefix"] = df["miRNA ID"].apply(get_prefix)
    logger.info(f"mirbase file after postprocess {df.shape}")
    df.to_csv(MIRBASE_FILE)
    logger.info("save mirbase file")



@click.group()
def cli():
    pass


# @click.command()
# @click.argument('min_version', type=float)
def mirbase_download(min_version: float):
    ftp = mirbase_connect()
    d = get_mirbase_versions(ftp, min_version)
    get_mature_file_from_dir_list(ftp, d)
    # read_fasta_files(min_version)
    # postprocess()

#
# cli.add_command(mirbase_download)
# cli.add_command(read_fasta_files)
# cli.add_command(postprocess)


if __name__ == '__main__':
    # cli()
    mirbase_download(17)

