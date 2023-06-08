from pathlib import Path

import click
from pandas import DataFrame

from utils.logger import logger
from utils.utilsfile import read_csv, to_csv
import pandas as pd


def read_blast_result_file(fin: Path) -> DataFrame:
    logger.info(f"read_blast_result_file {fin}")
    df: DataFrame = read_csv(Path(fin))
    df = df.astype({'s.start': 'Int32',
                    's.end': 'Int32'},
                   errors="ignore")

    region = fin.stem.split("_")[-1]
    df["region"] = region

    # take only the rows with valid results
    df.dropna(axis=0, how='any', subset=['sequence'], inplace=True)
    return df

@click.command()
@click.argument('directory', type=Path)
@click.argument('fname', type=str)
@click.argument('blast_prev_step_file', type=Path)
@click.argument('fout', type=Path)
def concat_blast_result(directory: Path, fname: str, blast_prev_step_file: Path, fout: Path):

    blast_result_list = [read_blast_result_file(f) for f in directory.glob(f"*{fname}_*.csv")]
    logger.info("Finish read the files. start to concatenate")
    blast_result_df = pd.concat(blast_result_list, axis=0, ignore_index=True)

    vc = blast_result_df["key"].value_counts()
    blast_result_df["region count"] = \
        blast_result_df.merge(vc, how="left", left_on="key", right_index=True)["key_y"]
    blast_result_inx = blast_result_df["key"].unique()

    all_interactions: DataFrame = read_csv(blast_prev_step_file)
    all_interactions["region"] = "None"
    all_interactions["region count"] = 0
    all_interactions.query("key not in @blast_result_inx", inplace=True)

    unite = pd.concat([blast_result_df, all_interactions], axis=0)
    unite.sort_values(by="key", ignore_index=True, inplace=True)
    unite.drop(columns=["start", "end"], inplace=True, errors="ignore")
    unite = unite.rename(columns={"s.start": "start",
                                  "s.end": "end"})
    to_csv(unite, fout)





@click.group()
def cli():
    pass


cli.add_command(concat_blast_result)



if __name__ == '__main__':
    cli()









