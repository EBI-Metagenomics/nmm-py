from ._alphabet import Alphabet, CAlphabet
from typing import Dict, Union

from ._ffi import ffi, lib


class CCodonTable:
    """
    Wrapper around the C implementation of a codon table.

    Parameters
    ----------
    imm_codont : ffi.CData
        Passing `None` will create a new codon table at the underlying library level using the
        `alphabet` argument.
    alphabet : Union[CAlphabet, None]
        Passing a `CAlphabet` will create a new codon table at the underlying library level.
    """

    def __init__(
        self, imm_codont: Union[ffi.CData, None], alphabet: Union[CAlphabet, None]
    ):
        if imm_codont is None:
            if alphabet is None:
                raise ValueError("`alphabet` is `None`")

            self.__cdata = lib.nmm_codont_create(alphabet.imm_abc)
        else:
            if alphabet is not None:
                raise ValueError("`alphabet` is not `None`")

            self.__cdata = imm_codont

        if self.__cdata == ffi.NULL:
            raise RuntimeError("`cdata` is NULL.")

    @property
    def nmm_codon(self) -> ffi.CData:
        return self.__cdata

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_codont_get_abc(self.__cdata)

    def set_lprob(self, seq: bytes, lprob: float) -> None:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        err: int = lib.nmm_codont_set_lprob(self.__cdata, create_codon(seq), lprob)
        if err != 0:
            s = seq.decode()
            raise ValueError(f"Could not set a probability for `{s}`.")

    def get_lprob(self, seq: bytes) -> float:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        return lib.nmm_codont_get_lprob(self.__cdata, create_codon(seq))

    def normalize(self) -> None:
        err: int = lib.nmm_codont_normalize(self.__cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self.__cdata != ffi.NULL:
            lib.nmm_codont_destroy(self.__cdata)


class CodonTable(CCodonTable):
    """
    Codon table.

    It is a table of probabilities for codons, in log-space.

    Parameters
    ----------
    alphabet : Alphabet.
    lprobs : Emission probabilities in log-space.
    """

    def __init__(self, alphabet: Alphabet, lprobs: Dict[bytes, float] = {}):
        super().__init__(None, alphabet)
        self._alphabet = alphabet
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


def create_codon(seq: bytes):
    codon = ffi.new("struct nmm_codon *")
    codon.a = seq[0:1]
    codon.b = seq[1:2]
    codon.c = seq[2:3]
    return codon
