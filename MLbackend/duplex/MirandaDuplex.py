import os
from tempfile import NamedTemporaryFile
from typing import List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from duplex.Duplex import Duplex, NoDuplexResult


def run_miranda(mirna: str, target: str) -> List:
    with NamedTemporaryFile(prefix="miranda_mrna_") as mrna_fasta_filename:
        with NamedTemporaryFile(prefix="miranda_mirna_") as mirna_fasta_filename:
            with NamedTemporaryFile(prefix="miranda_out_") as miranda_out:
                mRNA_record = SeqRecord(Seq(target), description="mRNA")
                miRNA_record = SeqRecord(Seq(mirna), description="miRNA")
                SeqIO.write(mRNA_record, mrna_fasta_filename.name, "fasta")
                SeqIO.write(miRNA_record, mirna_fasta_filename.name, "fasta")

                miranda_cmd = "miranda {mir} {mrna} -out {out} -en 10000 -sc 60 ".format(mir=mirna_fasta_filename.name,
                                                                                         mrna=mrna_fasta_filename.name,
                                                                                         out=miranda_out.name)
                os.system(miranda_cmd)

                with open(miranda_out.name, "r") as fle:
                    mr = fle.readlines()
    # for m in mr:
    #     print (m)
    return mr

def extract_seq (s):
        return s.split("'")[1].strip()[0:-2]

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]


class MirandaDuplex(Duplex):
    @classmethod
    def createDuplex(cls, mirna: str, target: str) -> Tuple[str, str, str, str]:
        mr = run_miranda(mirna, target)
        query_list = lines_that_contain("Query:", mr)
        info_list = lines_that_contain(">", mr)
        interaction_list = lines_that_contain("|", mr)
        ref_list = lines_that_contain("Ref:", mr)

        if len(query_list)<1:
            raise NoMirandaHits

        query = extract_seq(query_list[0].strip("\n"))
        ref = extract_seq(ref_list[0].strip("\n"))
        interaction = interaction_list[0][query_list[0].find(query):query_list[0].find(query)+len(query)]

        mrna_bulge = ""
        mrna_inter = ""
        mir_inter = ""
        mir_bulge = ""

        query = query[::-1].upper()
        ref = ref[::-1].upper()
        interaction = interaction[::-1]

        for i in range (len(interaction)):
            if interaction[i] == " " :
                mrna_inter += " "
                mir_inter += " "
                mrna_bulge += ref[i]
                mir_bulge += query[i]
            else :
                mrna_bulge += " "
                mrna_inter += ref[i]
                mir_inter += query[i]
                mir_bulge += " "
        mrna_bulge = mrna_bulge.replace("-", " ")
        mrna_inter = mrna_inter.replace("-", " ")
        mir_inter = mir_inter.replace("-", " ")
        mir_bulge = mir_bulge.replace("-", " ")

        return mrna_bulge, mrna_inter, mir_inter, mir_bulge


class NoMirandaHits(NoDuplexResult):
    pass


if __name__ == '__main__':
    mirna = 'UGAGGUAGUAGGUUGUAUAGUU'
    target = 'CAACAAGAACCAACACUACUGCCCAACCGUCAACGUCG'
    dp = MirandaDuplex.fromChimera(mirna, target)
    print(dp)
