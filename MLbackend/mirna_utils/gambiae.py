import pandas as pd
from Bio import SeqIO
from pandas import DataFrame

from consts.mirna_utils import GAMBIAE_FILE, MIRBASE_FILE
from utils.logger import logger
from utils.utilsfile import fasta_to_dataframe


def read_aga_fasta():
    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE, index_col=0)

    ver = "aga_paper"
    gambiae_mirna: DataFrame = fasta_to_dataframe(GAMBIAE_FILE, "aga-")
    gambiae_mirna.insert(0, "ver", ver)

    assert list(mirbase_df.columns) == list(gambiae_mirna.columns), \
        f"Can't concatenate the dataframe since they don't have the same columns name: " \
        f"{list(mirbase_df.columns)}" \
        f"{list(gambiae_mirna.columns)}"

    merge_df = pd.concat([mirbase_df, gambiae_mirna], ignore_index=True)
    merge_df.to_csv(MIRBASE_FILE)
    logger.info("update mirbase file")


if __name__ == '__main__':
    read_aga_fasta()

