from pathlib import Path

from consts.global_consts import ROOT_PATH, DATA_PATH

DATA_DIR = DATA_PATH / "mirna_utils/data"
FASTA_DIR = DATA_DIR / "fasta"
MIRBASE_URL = "mirbase.org"
VERSIONS_DIR = "pub/mirbase/"
MATURE_FILE = "mature.fa.gz"
MIRBASE_FILE: Path = DATA_DIR /"mirbase.csv"
GAMBIAE_FILE: Path = DATA_DIR / "aga_matureUpdate032616.fa"