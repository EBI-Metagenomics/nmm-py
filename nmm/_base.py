from ._alphabet import Alphabet
from typing import Dict

from ._ffi import ffi, lib


class CBase:
    def __init__(self, cdata: ffi.CData):
        self.__cdata = cdata

    @property
    def nmm_base(self) -> ffi.CData:
        return self.__cdata

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_base_get_abc(self.__cdata)

    def set_lprob(self, nucleotide: bytes, lprob: float) -> None:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        err: int = lib.nmm_base_set_lprob(self.__cdata, letter, lprob)
        if err != 0:
            nucl = nucleotide.decode()
            raise ValueError(f"Could not set a probability for `{nucl}`.")

    def get_lprob(self, nucleotide: bytes) -> float:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        return lib.nmm_base_get_lprob(self.__cdata, letter)

    def normalize(self) -> None:
        err: int = lib.nmm_base_normalize(self.__cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self.__cdata != ffi.NULL:
            lib.nmm_base_destroy(self.__cdata)


class Base(CBase):
    def __init__(self, alphabet: Alphabet, lprobs: Dict[bytes, float] = {}):
        self._alphabet = alphabet
        cdata = lib.nmm_base_create(alphabet.imm_abc)
        super().__init__(cdata)
        for letter, lprob in lprobs.items():
            self.set_lprob(letter, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet
