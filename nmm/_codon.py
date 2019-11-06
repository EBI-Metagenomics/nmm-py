from ._alphabet import Alphabet
from typing import Dict

from ._ffi import ffi, lib


class CCodon:
    def __init__(self, cdata: ffi.CData):
        self.__cdata = cdata

    @property
    def nmm_codon(self) -> ffi.CData:
        return self.__cdata

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_codon_get_abc(self.__cdata)

    def set_lprob(self, seq: bytes, lprob: float) -> None:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        a = seq[0:1]
        b = seq[1:2]
        c = seq[2:3]
        err: int = lib.nmm_codon_set_lprob(self.__cdata, a, b, c, lprob)
        if err != 0:
            s = seq.decode()
            raise ValueError(f"Could not set a probability for `{s}`.")

    def get_lprob(self, seq: bytes) -> float:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        return lib.nmm_codon_get_lprob(self.__cdata, seq[0:1], seq[1:2], seq[2:3])

    def normalize(self) -> None:
        err: int = lib.nmm_codon_normalize(self.__cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self.__cdata != ffi.NULL:
            lib.nmm_codon_destroy(self.__cdata)


class Codon(CCodon):
    def __init__(self, alphabet: Alphabet, lprobs: Dict[bytes, float] = {}):
        self._alphabet = alphabet
        cdata = lib.nmm_codon_create(self._alphabet.imm_abc)
        super().__init__(cdata)
        for seq, lprob in lprobs.items():
            self.set_lprob(seq, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet
