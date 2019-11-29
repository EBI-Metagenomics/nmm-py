from ._alphabet import Alphabet, CAlphabet
from typing import Dict, Optional

from ._ffi import ffi, lib


class CCodonTable:
    """
    Wrapper around the C implementation of a codon table.

    Parameters
    ----------
    nmm_codont : Optional[ffi.CData]
        Passing `None` will create a new codon table at the underlying library level using the
        `alphabet` argument.
    alphabet : Optional[CAlphabet]
        Passing a `CAlphabet` will create a new codon table at the underlying library level.
    """

    def __init__(
        self,
        nmm_codont: Optional[ffi.CData] = None,
        alphabet: Optional[CAlphabet] = None,
    ):
        self.__nmm_codont = ffi.NULL
        if nmm_codont is None:
            if alphabet is None:
                raise ValueError("`alphabet` is `None`")

            self.__nmm_codont = lib.nmm_codont_create(alphabet.imm_abc)
        else:
            if alphabet is not None:
                raise ValueError("`alphabet` is not `None`")

            self.__nmm_codont = nmm_codont

        if self.__nmm_codont == ffi.NULL:
            raise RuntimeError("`cdata` is NULL.")

    @property
    def nmm_codont(self) -> ffi.CData:
        return self.__nmm_codont

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_codont_get_abc(self.__nmm_codont)

    def set_lprob(self, seq: bytes, lprob: float) -> None:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        err: int = lib.nmm_codont_set_lprob(self.__nmm_codont, create_codon(seq), lprob)
        if err != 0:
            s = seq.decode()
            raise ValueError(f"Could not set a probability for `{s}`.")

    def get_lprob(self, seq: bytes) -> float:
        if len(seq) != 3:
            raise ValueError("Codon must have three letters.")

        return lib.nmm_codont_get_lprob(self.__nmm_codont, create_codon(seq))

    def normalize(self) -> None:
        err: int = lib.nmm_codont_normalize(self.__nmm_codont)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self.__nmm_codont != ffi.NULL:
            lib.nmm_codont_destroy(self.__nmm_codont)


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
        super().__init__(alphabet=alphabet)
        self._alphabet = alphabet
        for seq, lprob in lprobs.items():
            self.set_lprob(seq, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet


def create_codon(seq: bytes):
    codon = ffi.new("struct nmm_codon *")
    codon.a = seq[0:1]
    codon.b = seq[1:2]
    codon.c = seq[2:3]
    return codon
