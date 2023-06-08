import os
from tempfile import NamedTemporaryFile
from typing import List, Tuple

import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from duplex.Duplex import Duplex


def run_rnahybrid(mirna: str, target: str) -> List:
    cmd = "RNAhybrid -s 3utr_human {} {}".format(target, mirna)
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    rh = p.stdout.readlines()
    rh = [l.decode(encoding='unicode-escape') for l in rh]
    # for l in rh:
    #     print (l)
    return rh


class rnaHybrid(Duplex):
    @classmethod
    def createDuplex(cls, mirna: str, target: str) -> Tuple[str, str, str, str]:
        a = run_rnahybrid(mirna, target)
        b = [[a[l], a[l + 1], a[l + 2], a[l + 3]] for l in range(len(a)) if a[l].startswith('target 5')]
        end = b[0][0].find("3") - 1
        f = [t[10:end] for t in b[0]]
        f = [t[::-1] for t in f]

        mrna_bulge, mrna_inter, mir_inter, mir_bulge = f
        return mrna_bulge, mrna_inter, mir_inter, mir_bulge


    #     print("FFFFFFFFFFFFFFFFFF")
    #     print(f)
    #     tar = interaction_formatting(f[0], f[1])[1:]
    #     tar += "-" * max(0, (22 - len(tar)))
    #     # tar = tar[:22]
    #
    #     mir = interaction_formatting(f[3], f[2])[1:]
    #     mir += "-" * max(0, (22 - len(mir)))
    #     # mir = mir[:22]
    #     print(mir)
    #     print(tar)
    #     exit(3)
    # #     return tar, mir
    # #
    # # return mrna_bulge, mrna_inter, mir_inter, mir_bulge

#
# class NoMirandaHits(Exception):
#     pass


if __name__ == '__main__':
    mirna = 'UGAGGUAGUAGGUUGUAUAGUU'
    target = 'CAACAAGAACCAACACUACUGCCCAACCGUCAACGUCG'
    dp = MirandaDuplex.fromChimera(mirna, target)
    print(dp)
