from functools import partial
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Tuple

import numpy as np
import pandas as DataFrame
import pandas as pd
from Bio import SeqIO


import os
os.environ['PATH'] = '/home/efrco/.conda/envs/my_env/bin:/storage/modules/packages/anaconda3/bin:/storage/modules/bin:/storage/modules/packages/anaconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/storage/modules/packages/matlab/R2019B/bin:/home/efrco/.local/bin:/home/efrco/bin'

from Bio.Blast.Applications import NcbiblastnCommandline

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pandas import Series

from consts.biomart import BIOMART_BLAST_PATH, BIOMART_DATA_PATH
from consts.global_consts import ROOT_PATH
from logger import logger

from utilsfile import call_wrapper, read_csv, get_wrapper, to_csv

from consts.global_consts import MINIMAL_LENGTH_TO_BLAST, MINIMAL_BLAST_IDENTITY, MINIMAL_BLAST_COVERAGE
from consts.pipeline_steps import REGION_PATH


def run_blastn(seq: str, db_title: str) -> Series:
    def blast_coverage(start: int, end: int, query: str) -> float:
        return (end - start + 1.0) * 100 / len(query)

    RETURN_COL = ["Gene_ID", "sequence", "identity", "coverage", "s.start", "s.end"]
    try:
        if len(seq) < MINIMAL_LENGTH_TO_BLAST:
            return pd.Series(index=RETURN_COL)
    except TypeError:
        return pd.Series(index=RETURN_COL)


    with NamedTemporaryFile(prefix="blast") as blast_out_file:
        with NamedTemporaryFile(prefix="blast") as seq_to_find_file:
            record = SeqRecord(Seq(seq), description="seq_to_find")

            SeqIO.write(record, seq_to_find_file.name, "fasta")

            cline = NcbiblastnCommandline(query=str(seq_to_find_file.name), db=db_title, evalue=1,
                                          strand="plus", task="blastn-short",
                                          out=str(blast_out_file.name), outfmt=6)
            print(cline)

            call_wrapper(cmd=str(cline), cwd=BIOMART_BLAST_PATH)

        #Parse the output file
        colnames = ['query acc.ver', 'subject acc.ver', '%identity', 'alignment length', 'mismatches',
                    'gap opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit score']
        result = pd.read_csv(blast_out_file.name, sep='\t', names=colnames, header=None)
        result.rename(columns={'%identity': "identity",
                               'alignment length': 'alignment_length',
                               "gap opens" : "gap_opens"}, inplace=True)


        try:
            result["coverage"] = result.apply(func=get_wrapper(blast_coverage,
                                                               's.start', 's.end',
                                                               query=seq),
                                              axis=1)
        except ValueError:
            # Empty result dataframe
            assert result.shape[0] == 0, "Wrong exception. have to check"
            return pd.Series(index=RETURN_COL)


        # Consider the full match rows only
        result.query("identity >= @MINIMAL_BLAST_IDENTITY and "
                     "coverage >= @MINIMAL_BLAST_COVERAGE and "
                     "gap_opens == 0", inplace=True)
        result.reset_index(inplace=True)

        if result.shape[0] == 0:
            return pd.Series(index=RETURN_COL)

        # get the full sequence
        transcripts: DataFrame = pd.read_csv(BIOMART_DATA_PATH / f"{db_title}.csv")
        result = result.merge(transcripts, how="left", left_on="subject acc.ver", right_on="ID")

        # choose the row with longest utr and add the full mrna
        ###############################
        best = result.iloc[result["sequence length"].idxmax()]
        best.rename({'ID': "Gene_ID"}, inplace=True)


        return best[RETURN_COL]
        # return str(best["sequence"].iloc[0]), best["s.start"].iloc[0], best["s.end"].iloc[0]


def get_blast_result(seq: str, db_title: str) -> Tuple[str, str]:
    seq = run_blastn(seq, db_title)
    region = "" if np.isnan else db_title
    return region, seq


def blast_file(fin: Path, fout: Path, db_title: str):

    logger.info(f"blast file {fin} against {db_title}")
    in_df: DataFrame = read_csv(fin)
    blastn_df: DataFrame = in_df.apply(func=get_wrapper(run_blastn,
                                                        "site", db_title=db_title),
                                       axis=1)
    result = pd.concat([in_df, blastn_df], axis=1)

    # in_df["blast region"] = in_df["blast sequence"].apply(lambda x: "" if np.isnan(x) else db_title)

    to_csv(result, fout)


def df_contains(substr: str, df: DataFrame) -> str:
    """
    find the the longest match row in dataframe
    """

    df = df[df["sequence"].str.contains(substr, regex=False)]
    if len(df) == 0:
        return ""
    df.sort_values(by="sequence length", ascending=False, inplace=True)
    return str(df.iloc[0]["sequence"])


def fast_blast_file(fin: Path, fout: Path, db_title: str):
    logger.info(f"fast blast file {fin} against {db_title}")
    in_df: DataFrame = read_csv(fin)
    seq_file = BIOMART_DATA_PATH / "human_3utr.csv"
    df_contains_db_title = partial(df_contains, df=pd.read_csv(seq_file))

    in_df["blast sequence"] = in_df.apply(func=get_wrapper(df_contains_db_title,
                                                           "full_mrna"),
                                          axis=1)
    to_csv(in_df, fout)


# fin = ROOT_PATH / "generate_interactions/new.csv"
# fout = ROOT_PATH / "generate_interactions/blast.csv"
# fast_blast_file(fin, fout, "human_3utr.csv")
#


def run_blastn(seq: str, db_title: str) -> Series:
    def blast_coverage(start: int, end: int, query: str) -> float:
        return (end - start + 1.0) * 100 / len(query)

    RETURN_COL = ["mRNA_name", "full_mrna", "identity", "coverage", "s.start", "s.end"]
    try:
        if len(seq) < MINIMAL_LENGTH_TO_BLAST:
            return pd.Series(index=RETURN_COL)
    except TypeError:
        return pd.Series(index=RETURN_COL)


    with NamedTemporaryFile(prefix="blast") as blast_out_file:
        with NamedTemporaryFile(prefix="blast") as seq_to_find_file:
            record = SeqRecord(Seq(seq), description="seq_to_find")

            SeqIO.write(record, seq_to_find_file.name, "fasta")
            print(seq_to_find_file.name)
            cline = NcbiblastnCommandline(query=str(seq_to_find_file.name), db=db_title, evalue=1,
                                          strand="plus", task="blastn-short",
                                          out=str(blast_out_file.name), outfmt=6)

            call_wrapper(cmd=str(cline), cwd=BIOMART_BLAST_PATH)

        #Parse the output file
        colnames = ['query acc.ver', 'subject acc.ver', '%identity', 'alignment length', 'mismatches',
                    'gap opens', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit score']
        result = pd.read_csv(blast_out_file.name, sep='\t', names=colnames, header=None)
        result.rename(columns={'%identity': "identity",
                               'alignment length': 'alignment_length',
                               "gap opens": "gap_opens"}, inplace=True)


        try:
            result["coverage"] = result.apply(func=get_wrapper(blast_coverage,
                                                               's.start', 's.end',
                                                               query=seq), axis=1)
        except ValueError:
            # Empty result dataframe
            assert result.shape[0] == 0, "Wrong exception. have to check"
            return pd.Series(index=RETURN_COL)


        # Consider the full match rows only
        result.query("identity >= @MINIMAL_BLAST_IDENTITY and "
                     "coverage >= @MINIMAL_BLAST_COVERAGE and "
                     "gap_opens == 0", inplace=True)
        result.reset_index(inplace=True)

        if result.shape[0] == 0:
            return pd.Series(index=RETURN_COL)

        # get the full sequence
        #transcripts: DataFrame = pd.read_csv(BIOMART_DATA_PATH / f"{db_title}.csv")
        transcripts: DataFrame = pd.read_csv(ROOT_PATH / "generate_interactions/new_2.csv")

        result = result.merge(transcripts, how="left", left_on="subject acc.ver", right_on="Gene_ID")

        # choose the row with longest utr and add the full mrna
        ###############################
        #best = result.iloc[result["sequence length"].idxmax()]
        #best.rename({'ID': "Gene_ID"}, inplace=True)


        return result
        # return str(best["sequence"].iloc[0]), best["s.start"].iloc[0], best["s.end"].iloc[0]

path= BIOMART_BLAST_PATH / "human_3utr"
r=run_blastn("AAGAAAGATACTCATTTATAGTTACGTTCATTTCAGGTTAAACATGAAAGAAGCCTGGTTACTGATTTGTATAAAATGTACTCTTAAAGTATAAAATATAAGGTAAGGTAAATTTCATGCATCTTTTTATGAAGACCACCTATTTTATATTTCAAATTAAATAATTTTAAAGTTGCTGGCCTAATGAGCAA"
              "TGTTCTCAATTTTCGTTTTCATTTTGCTGTATTGAGACCTATAAATAAATGTATATTTTTTTTTGCATAAAGTA",path)
print(r)