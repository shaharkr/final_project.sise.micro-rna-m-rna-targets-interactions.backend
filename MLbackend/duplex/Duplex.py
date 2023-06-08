from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Tuple

from consts.features import GU, MM
from duplex.utils import mix_inter_bulge_seq


class DuplexException(Exception):
    pass


class DuplexRepresentationException(DuplexException):
    """Both interaction and bulge have nt value. one of them should be empty (space)."""
    pass


class NoDuplexResult(DuplexException):
    pass


def padding(s: str, max_len: int) -> str:
    return s + " " * (max_len - len(s))



# class Duplex(ABC):
class Duplex():
    def __init__(self, mrna_bulge: str, mrna_inter: str, mir_inter: str, mir_bulge:str) -> None:
        max_len = max(map(len, [mrna_bulge, mrna_inter, mir_inter, mir_bulge]))
        self._mrna_bulge = padding(mrna_bulge, max_len)
        self._mrna_inter = padding(mrna_inter, max_len)
        self._mir_inter = padding(mir_inter, max_len)
        self._mir_bulge = padding(mir_bulge, max_len)

        self.__verify()

    @classmethod
    def fromStrings(cls, mrna_bulge: str, mrna_inter: str, mir_inter: str, mir_bulge: str) -> Duplex:
        return cls(mrna_bulge, mrna_inter, mir_inter, mir_bulge)

    @classmethod
    def fromChimera(cls, mirna: str, target: str) -> Duplex:
        try:
            mrna_bulge, mrna_inter, mir_inter, mir_bulge = cls.createDuplex(mirna, target)

            reconstruct_mrna = mix_inter_bulge_seq(mrna_bulge, mrna_inter)
            reconstruct_mirna = mix_inter_bulge_seq(mir_inter, mir_bulge)

            assert reconstruct_mrna[::-1].replace("#", "") in target, \
                f"""target:     {target}
    reconstruct: {reconstruct_mrna[::-1]}
    mrna_bulge: {mrna_bulge}
    mrna_inter: {mrna_inter}"""

            assert reconstruct_mirna in mirna.replace("#", ""), f"""mirna:    {mirna}
    reconstruct: {reconstruct_mirna}
    mir_inter:   {mir_inter}
    mir_bulge:   {mir_bulge}"""

        except NoDuplexResult:
            mrna_bulge, mrna_inter, mir_inter, mir_bulge = "", "", "", ""

        return cls(mrna_bulge, mrna_inter, mir_inter, mir_bulge)

    @classmethod
    @abstractmethod
    def createDuplex(cls, mirna: str, target: str) -> Tuple[str, str, str, str]:
        pass


    def serialize(self) -> Tuple[str, str, str, str]:
        return self._mrna_bulge, self._mrna_inter, self._mir_inter, self._mir_bulge

    @property
    def valid(self) -> bool:
        return self.mrna_bulge != "" and self.mrna_inter != "" and self.mir_inter != "" and self.mir_bulge != ""

    @property
    def mrna_bulge(self):
        return self._mrna_bulge

    @property
    def mrna_inter(self):
        return self._mrna_inter

    @property
    def mir_inter(self):
        return self._mir_inter

    @property
    def mir_bulge(self):
        return self._mir_bulge

    def __verify(self) -> None:
        def check_pairs(p1, p2):
            for i in range(min(len(p1), len(p2))):
                if p1[i] != ' ' and p2[i] != ' ':
                    raise DuplexRepresentationException

        check_pairs(self._mrna_bulge, self._mrna_inter)
        check_pairs(self._mir_bulge, self._mir_inter)

    def tostring(self) -> str:
        classstr = ""
        classstr = classstr + "target_bulge:       {}\n".format(self._mrna_bulge)
        classstr = classstr + "target_interaction: {}\n".format(self._mrna_inter)
        classstr = classstr + "mirna_interaction:  {}\n".format(self._mir_inter)
        classstr = classstr + "mirna_bulge:        {}\n".format(self._mir_bulge)
        classstr +="\n"
        # classstr = classstr + "site ({}):          {}\n".format(len(self.site), self.site)

        #classstr = classstr + "mrna_bulges_count: {} \nmir_bulges_count:  {}\n".format(self.mrna_bulges_count,self.mir_bulges_count)

        return classstr

    def __str__(self):
        return self.tostring()

    def interaction_iterator(self):
        for i in range(len(self._mir_inter)):
            if self._mir_inter[i] != ' ':
                yield i, self._mrna_inter[i] + self._mir_inter[i]

    def mir_iterator (self):
        #make sure the iterator won't access out of range of the mir variables
        mir_len = max(len(self._mir_inter), len(self._mir_bulge))
        mir_inter = self._mir_inter + " "*mir_len
        mir_bulge = self._mir_bulge + " "*mir_len

        for i in range(mir_len):
            if mir_inter[i]!=' ' or mir_bulge[i]!=' ':
                yield i, mix_inter_bulge_seq(mir_inter[i], mir_bulge[i])

    def pair_iterator(self):
        mir_len = max(len(self._mir_inter), len(self._mir_bulge))
        mir_inter = self._mir_inter + " " * (mir_len - len(self._mir_inter))
        mir_bulge = self._mir_bulge + " " * (mir_len - len(self._mir_bulge))

        for i in range(mir_len):
            mir = mix_inter_bulge_seq(mir_inter[i], mir_bulge[i])
            try:
                m_i = self._mrna_inter[i]
            except IndexError:
                m_i = " "
            try:
                m_b = self._mrna_bulge[i]
            except IndexError:
                m_b = " "

            mrna = mix_inter_bulge_seq(m_i, m_b)

            yield mrna+mir


    # def mir_pairing_iterator (self):
    #     for i, mir in self.mir_iterator():
    #         try:
    #             m_i = self._mrna_inter[i]
    #         except IndexError:
    #             m_i = " "
    #         try:
    #             m_b = self._mrna_bulge[i]
    #         except IndexError:
    #             m_b = " "
    #
    #         mrna = mix_inter_bulge_seq(m_i, m_b)
    #
    #         yield mrna+mir

    def extract_seed(self, start: int, end: int) -> Duplex:
        # notice: if there a bulge in the mrna, it will skip it. the function return 8 consecutive mirna nt.
        mir_i = self.mir_iterator()

        for i in range(start - 1):
            next(mir_i)
        s, _ = next(mir_i)

        e = None
        for i in range(end - start):
            e, _ = next(mir_i)
        e += 1

        return self.fromStrings(self._mrna_bulge[s:e],
                                self._mrna_inter[s:e],
                                self._mir_inter[s:e],
                                self._mir_bulge[s:e])

    @property
    def seed(self) -> Duplex:
        return self.extract_seed(start=1, end=8)

    @property
    def site(self) -> str:
        def first_char(s1, s2):
            i = 0
            while s1[i] == " " and s2[i] == " ":
                i+=1
            return i

        # make sure the iterator won't access out of range of the mir variables
        mrna_site = ""
        mir_len = max(len(self._mir_inter), len(self._mir_bulge))

        mrna_inter = self._mrna_inter + " "*mir_len
        mrna_bulge = self._mrna_bulge + " "*mir_len

        mir_start = first_char(self._mir_inter, self._mir_bulge)
        for i in range(mir_start, mir_len):
            mrna_site += mix_inter_bulge_seq(mrna_inter[i], mrna_bulge[i])

        mrna_site = mrna_site.replace(" ","")


        # print("DUPLEX:", len(mrna_site))
        # print("DUPLEX:", mrna_site)
        return mrna_site

    @property
    def mir_bulge_count(self) -> int:
        mir_bulge  = self._mir_bulge + " "*len(self._mrna_bulge)
        mrna_bulge = self._mrna_bulge + " "*len(self._mir_bulge)

        a = ""
        for i in range(len(mir_bulge)):
            a += "-" if mir_bulge[i] != " " and mrna_bulge[i] == " " else " "
        return len(a.split())

    @property
    def mrna_bulge_count(self) -> int:
        mir_bulge = self._mir_bulge + " " * len(self._mrna_bulge)
        mrna_bulge = self._mrna_bulge + " " * len(self._mir_bulge)

        a = ""
        for i in range(len(mir_bulge)):
            a += "-" if mir_bulge[i] == " " and mrna_bulge[i] != " " else " "
        return len(a.split())

    @property
    def interaction_count(self) -> int:
        return sum(1 for _ in self.interaction_iterator())


    @property
    def canonical_seed(self):
        """canonical seed: exact W–C pairing of 2–7 or 3–8 nts of the miRNA"""
        seed: Duplex = self.extract_seed(1, 8)
        seed_2_7: Duplex = self.extract_seed(2, 7)
        seed_3_8: Duplex = self.extract_seed(3, 8)

        c2_7 = len([pair for _, pair in seed_2_7.interaction_iterator() if pair not in GU])
        c3_8 = len([pair for _, pair in seed_3_8.interaction_iterator() if pair not in GU])
        return (c2_7 == 6) or (c3_8 == 6)

    @property
    def noncanonical_seed(self):
        """non-canonical seed: pairing at positions 2–7 or 3–8,
        allowing G-U pairs and up to one bulged or mis-matched nucleotide"""

        seed: Duplex = self.extract_seed(1, 8)
        seed_2_7: Duplex = self.extract_seed(2, 7)
        seed_3_8: Duplex = self.extract_seed(3, 8)

        c2_7 = len([pair for _, pair in seed_2_7.interaction_iterator()])
        c3_8 = len([pair for _, pair in seed_3_8.interaction_iterator()])

        c2_7_MM = len([pair for pair in seed_2_7.pair_iterator() if pair in MM])
        c3_8_MM = len([pair for pair in seed_3_8.pair_iterator() if pair in MM])
        c2_7_B = len([pair for pair in seed_2_7.pair_iterator() if len(pair.strip()) == 1])
        c3_8_B = len([pair for pair in seed_3_8.pair_iterator() if len(pair.strip()) == 1])

        bulge_mismatch_2_7 = c2_7_MM + c2_7_B
        bulge_mismatch_3_8 = c3_8_MM + c3_8_B

        r2_7 = (c2_7 >= 5) and (bulge_mismatch_2_7 <= 1)
        r3_8 = (c3_8 >= 5) and (bulge_mismatch_3_8 <= 1)

        return r2_7 or r3_8

    @property
    def site_non_match_tail(self) -> str:

        mir_len = max(len(self._mir_inter), len(self._mir_bulge))
        mrna_bulge = self._mrna_bulge + " "*mir_len
        ca = self._mir_bulge[0]
        count_not_match = 0
        # check left side
        for i in range(0, len(self._mir_bulge)):
            if self._mir_bulge[i] != " ":
                if mrna_bulge[i] == " ":
                    count_not_match = count_not_match + 1
                else:
                    break
            else:
                break

        # check right side
        for i in range(len(self._mir_bulge) - 1,0,-1):
            if self._mir_bulge[i] != " ":
                if mrna_bulge[i] == " ":
                    count_not_match = count_not_match + 1
                else:
                    break
            else:
                break

        return count_not_match