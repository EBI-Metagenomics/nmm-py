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

        err: int = lib.nmm_codon_set_lprob(self.__cdata, ccodon_code(seq), lprob)
        if err != 0:
            s = seq.decode()
            raise ValueError(f"Could not set a probability for `{s}`.")

    def get_lprob(self, seq: bytes) -> float:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        return lib.nmm_codon_get_lprob(self.__cdata, ccodon_code(seq))

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


# class CCodonCode:
#     """
#     Wrapper around the C implementation of codon code.

#     Parameters
#     ----------
#     cdata : `struct imm_abc *`.
#     """

#     def __init__(self, cdata: ffi.CData):
#         self.__cdata = cdata

# def decode(self, seq: bytes) -> DecodedCodon:
#     ccode = ffi.new("struct nmm_ccode *")
#     lprob: float = lib.nmm_frame_state_decode(self._cdata, seq, len(seq), ccode)
#     return DecodedCodon(lprob, ccode.a + ccode.b + ccode.c)


def ccodon_code(seq: bytes):
    ccode = ffi.new("struct nmm_ccode *")
    ccode.a = seq[0:1]
    ccode.b = seq[1:2]
    ccode.c = seq[2:3]
    return ccode
