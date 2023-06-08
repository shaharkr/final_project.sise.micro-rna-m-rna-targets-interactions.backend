from functools import lru_cache
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from utilsfile import get_subsequence_by_coordinates


def read_genome(self):
    gen = SeqIO.parse(self.genome_file, "fasta")
    self.genome_list = list(gen)


@lru_cache(maxsize=None)
def get_chr(name: str, directory: Path):
    full_path: Path = directory / f"{name}.fa"
    seq_record = SeqIO.parse(full_path, "fasta")
    return next(seq_record).seq


def extract_seq_from_chromosome(chr_name: str, start: int, stop: int, strand: str, directory: Path,  extra_chars: int = 0) -> str:
    chr_seq = get_chr(chr_name, directory)
    seq = get_subsequence_by_coordinates(str(chr_seq), start, stop, strand, extra_chars)

    return seq

