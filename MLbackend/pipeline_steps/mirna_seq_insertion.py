from pathlib import Path
import pandas as pd
from pandas import DataFrame
import click

from consts.mirna_utils import MIRBASE_FILE
from consts.pipeline_steps import READ_PATH, MIRNA_SEQ_PATH
from utils.logger import logger
from utils.utilsfile import get_wrapper, read_csv, to_csv


@click.command()
@click.argument('fname', type=str)
def mirna_seq_insertion(fname: str):
    logger.info(f"Insert mirna sequence to {fname}")

    fin_full_path = READ_PATH / fname
    fout_full_path = MIRNA_SEQ_PATH / fname

    df: DataFrame = read_csv(fin_full_path)
    df.drop(columns=["miRNA sequence"], inplace=True, errors='ignore') #drop the col we want to add via the join
    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE,
                                        usecols=["miRNA ID", "miRNA sequence"])
    join_df = df.merge(mirbase_df, how="left", left_on="miRNA ID", right_on="miRNA ID")
    to_csv(join_df, fout_full_path)
    logger.info(f"Finish the mirna sequence insertion to {fname}")


def qclash_mirna_func(mirna_hairpin_id: str, mrina_sub_seq: str, mirbase_hsa: DataFrame) -> str:
    sep = "" if ("-5p" in mirna_hairpin_id or "-3p" in mirna_hairpin_id) else "-"
    candidates = mirbase_hsa[mirbase_hsa['miRNA ID'].str.contains(f"{mirna_hairpin_id}{sep}")]['miRNA sequence']
    candidates = candidates[candidates.str.contains(mrina_sub_seq)]
    try:
        return candidates.iloc[0]
    except IndexError:
        return "Error: No mirna"



@click.command()
@click.argument('fname', type=str)
def qclash_melanoma_mirna_seq_insertion(fname: str):
    logger.info(f"Insert mirna sequence to {fname}")

    fin_full_path = READ_PATH / fname
    fout_full_path = MIRNA_SEQ_PATH / fname

    df: DataFrame = read_csv(fin_full_path)
    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE,
                                        usecols=["miRNA ID", "miRNA sequence", "prefix"])

    hsa = mirbase_df.query("prefix == 'hsa'")
    df["miRNA sequence"] = df.apply(func=get_wrapper(qclash_mirna_func,
                                                     'miRNA ID', 'mirna_seq_tmp', mirbase_hsa=hsa),
                                    axis=1, result_type="expand")

    to_csv(df, fout_full_path)
    logger.info(f"Finish the mirna sequence insertion to {fname}")



@click.command()
@click.argument('fname', type=str)
def mirnaid_fix(fname: str):
    d = {"mouse": "mmu",
         "human": "hsa",
         "elegans": "cel",
         "cattle": "bta",
         "fly": "aga"
         }
    prefix = None
    for k,v in d.items():
        if k in fname:
            prefix = v
    if prefix is None:
        raise Exception("unrecognized mirbase prefix")

    mirbase_df: DataFrame = pd.read_csv(MIRBASE_FILE).query("prefix==@prefix")
    mirbase_df.sort_values(by="version", ascending=False, inplace=True)
    mirbase_df.drop_duplicates("miRNA sequence", keep="first", inplace=True)

    fin_full_path = READ_PATH / fname
    fout_full_path = MIRNA_SEQ_PATH / fname

    d: DataFrame = read_csv(fin_full_path)


    join_df = d.merge(mirbase_df, how="left", left_on="miRNA sequence", right_on="miRNA sequence")
    d['miRNA ID'] = join_df['miRNA ID_y']
    to_csv(d, fout_full_path)


@click.group()
def cli():
    pass


cli.add_command(mirna_seq_insertion)
cli.add_command(qclash_melanoma_mirna_seq_insertion)
cli.add_command(mirnaid_fix)



if __name__ == '__main__':
    cli()
