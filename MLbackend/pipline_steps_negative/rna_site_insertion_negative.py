from pathlib import Path
from pandas import DataFrame, Series
import click

from consts.global_consts import HUMAN_SITE_EXTENDED_LEN, SITE_EXTRA_CHARS
from consts.pipeline_steps import ROOT_PATH
# from utils.logger import logger
from utils.utilsfile import get_subsequence_by_coordinates, get_wrapper, read_csv, to_csv


# @click.command()
# @click.argument('fin', type=str)
# @click.argument('fout', type=str)
def insert_site_by_coordinates(fin: str, fout:str):
    # logger.info(f"Insert site to {fin}")
    df: DataFrame = read_csv(Path(fin))
    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "mRNA sequence", "chimera_start", "chimera_end",
                                           extra_chars=SITE_EXTRA_CHARS),
                          axis=1)

    to_csv(df, Path(fout))
    # logger.info(f"finish the site sequence insertion to {fin}")



def get_site_from_extended_site(fin: DataFrame):
    def calc_chimera_start(seq: str, subseq: str) -> int:
        subseq = subseq.replace("#", "")
        try:
            if seq.find(subseq) == -1:
                print(seq)
                print(subseq)
                return -1
            return seq.find(subseq) + 1
        except AttributeError:
            return -1

    def calc_chimera_end(start: int, site: str) -> int:
        site = site.replace("#", "")
        if start == -1:
            return -1
        return start + len(site) - 1

    df = fin
    df["start"] = df.apply(func=get_wrapper(calc_chimera_start,'full_mrna', 'site'), axis=1)

    df["end"] = df.apply(func=get_wrapper(calc_chimera_end,'start', 'site'),axis=1)

    # df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
    #                                        "full_mrna", "start", "end", extra_chars=SITE_EXTRA_CHARS),axis=1)

    df["site"] = df.apply(func=get_wrapper(get_subsequence_by_coordinates,
                                           "full_mrna", "start", "end"), axis=1)
    return df
    




if __name__ == '__main__':
    pass
    # fin = ROOT_PATH / "generate_interactions/duplex_sample.csv"
    # fout = ROOT_PATH / "generate_interactions/duplex_rna_site.csv"
    # get_site_from_extended_site(fin, fout)

