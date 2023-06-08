from pathlib import Path
from typing import List


import click
import pandas as pd
from pandas import DataFrame, Series

from consts.biomart import BIOMART_DATA_PATH
from consts.pipeline_steps import GAMBIAE_INFORMATION_COLUMN_NAMES, GAMBIAE_INFORMATION_DTYPE, \
    GAMBIAE_INFORMATION_FILENAME, \
    GAMBIAE_INFORMATION_USECOLS, REGION_PATH, SITE_PATH
from utils.logger import logger
from utils.utilsfile import concatenate_biomart_df, get_wrapper, read_csv, to_csv


def add_gambiae_region_information(fin: Path) -> DataFrame:
    logger.info(f"enter to add_gambiae_region_information")
    gambiae_region_information: DataFrame = pd.read_csv(GAMBIAE_INFORMATION_FILENAME,
                                                        names=GAMBIAE_INFORMATION_COLUMN_NAMES,
                                                        usecols=GAMBIAE_INFORMATION_USECOLS,
                                                        dtype=GAMBIAE_INFORMATION_DTYPE,
                                                        delimiter="\t")

    in_df: DataFrame = read_csv(fin)
    join_df = in_df.merge(gambiae_region_information, how="left",
                          left_on="mRNA ID",  right_on="TRANSCRIPT_ID")
    logger.info(f"Finish enter the gambiae region information")
    return join_df


def get_region_ranges(len_5utr: int, len_cds: int, len_3utr: int) ->List:
    start = 0
    end = 0
    lens = [len_5utr, len_cds, len_3utr]
    range_list = []
    for i in range(3):
        end += lens[i]
        range_list.append(range(start, end))
        start += lens[i]
    return range_list

def check_trascript_length():


# assert len(transcript) == (len_5utr + len_cds + len_3utr), \
#     f"transcript and region length error: {len(transcript)} != {(len_5utr + len_cds + len_3utr)}"
# assert chimera_start in range(0, len(transcript)), f"chimera_start is not in proper range. " \
#                                                    f"chimera_start={chimera_start}, range=0..{len(transcript)}"
# assert chimera_end in range(0, len(transcript)), f"chimera_end is not in proper range. " \
#                                                    f"chimera_end={chimera_end}, range=0..{len(transcript)}"
    return




def find_gambiae_region(len_5utr: int, len_cds: int, len_3utr: int,
                        chimera_start: int, chimera_end: int) -> str:
    chimera_start -= 1  # python is zerobase
    chimera_end -= 1    # python is zerobase

    range_list = get_region_ranges(len_5utr, len_cds, len_3utr)

    # check the belonging of the chimera to each region
    result = []
    for region in range(3):
        result.append(
            (chimera_start in range_list[region]) or
            (chimera_end in range_list[region]))

    # if chimera start is within the 5utr and chimera end is within the 3utr, than the cds have to be true as well
    result[1] = result[1] or (result[0] and result[2])
    result: Series = pd.Series(data=result, index=['utr5', 'cds', 'utr3'])
    return "+".join(result.index[result])


def insert_gambiae_region(df) ->DataFrame:
    logger.info(f"enter to insert_gambiae_region")
    df["region"] = df.apply(func=get_wrapper(find_gambiae_region,
                                             "LEN_5UTR", "LEN_CDS", "LEN_3UTR",
                                             "chimera_start", "chimera_end"),
                            axis=1)
    return df


def find_gambiae_region_sequence(transcript: str, region: str, len_5utr: int, len_cds: int, len_3utr)\
        -> str:

    region_key = {'utr5': 0,
                  'cds' : 1,
                  'utr3' :2}
    if region not in region_key:
        return "None"

    range_list = get_region_ranges(len_5utr, len_cds, len_3utr)
    current_range = range_list[region_key[region]]
    return transcript[current_range.start:current_range.stop]


def insert_gambiae_region_sequence(df) ->DataFrame:
    logger.info(f"enter to insert_gambiae_region_sequence")
    df["region_sequence"] = df.apply(
        func=get_wrapper(find_gambiae_region_sequence,
                         "mRNA sequence", "region", "LEN_5UTR", "LEN_CDS", "LEN_3UTR"),
        axis=1)
    return df


def human_mapping_merge_by_name(fin: Path, fout: Path):
    def verify_sequence(seq: str, subseq: str) -> bool:
        try:
            return seq.find(subseq) != -1
        except AttributeError:
            return False

    in_df: DataFrame = read_csv(fin)
    in_df["join_key"] = in_df["mRNA ID"].apply(lambda x: "|".join(x.split("_")[0:2]))
    mRNA_df = concatenate_biomart_df("human")

    in_df = in_df.merge(mRNA_df, how="left",
                left_on=["region", "join_key"],
                right_on=["region", "ID"])

    in_df = in_df.rename(columns={"sequence": "region sequence"})
    in_df = in_df[['key', 'paper name', 'miRNA ID', 'miRNA sequence', 'mRNA ID',
                   'mRNA_seq_extended', 'region', 'region_sequence', 'mRNA_start', 'mRNA_end_extended']]


    in_df["join_ok"] = in_df.apply(func=get_wrapper(verify_sequence, 'region sequence', 'mRNA_seq_extended'), axis=1)

    to_csv(in_df, fout)





@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def gambiae_run(fin: str, fout:str):
    df: DataFrame = add_gambiae_region_information(Path(fin))
    df = insert_gambiae_region(df)
    df = insert_gambiae_region_sequence(df)
    df["start"] = df.apply(lambda row: row['chimera_start'] - row['mRNA sequence'].find(row['region_sequence']), axis=1)
    df["end"] = df.apply(lambda row: row['chimera_end'] - row['mRNA sequence'].find(row['region_sequence']), axis=1)
    df.rename(columns={"TRANSCRIPT_ID": ",Gene_ID"}, inplace=True)
    cols = [c for c in df.columns if c not in GAMBIAE_INFORMATION_USECOLS]

    to_csv(df[cols], Path(fout))


@click.command()
@click.argument('fin', type=str)
@click.argument('fout', type=str)
def human_mapping_run(fin: str, fout: str):
    human_mapping_merge_by_name(Path(fin), Path(fout))


    # df: DataFrame = add_gambiae_region_information(Path(fin))
    # df = insert_gambiae_region(df)
    # df = insert_gambiae_region_sequence(df)
    #
    # cols = [c for c in df.columns if c not in GAMBIAE_INFORMATION_USECOLS]
    # to_csv(df[cols], Path(fout))


@click.group()
def cli():
    pass


cli.add_command(gambiae_run)
cli.add_command(human_mapping_run)



if __name__ == '__main__':
    cli()
