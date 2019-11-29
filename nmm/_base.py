from ._alphabet import Alphabet, CAlphabet
from typing import Dict, Optional

from ._ffi import ffi, lib


class CBaseTable:
    """
    Wrapper around the C implementation of a base table.

    Parameters
    ----------
    nmm_baset : Optional[ffi.CData]
        Passing `None` will create a new base at the underlying library level using the `alphabet`
        argument.
    alphabet : Optional[CAlphabet]
        Passing a `CAlphabet` will create a new base at the underlying library level.
    """

    def __init__(
        self,
        nmm_baset: Optional[ffi.CData] = None,
        alphabet: Optional[CAlphabet] = None,
    ):
        if nmm_baset is None:
            if alphabet is None:
                raise ValueError("`alphabet` is `None`")

            self.__nmm_baset = lib.nmm_baset_create(alphabet.imm_abc)
        else:
            if alphabet is not None:
                raise ValueError("`alphabet` is not `None`")

            self.__nmm_baset = nmm_baset

        if self.__nmm_baset == ffi.NULL:
            raise RuntimeError("`cdata` is NULL.")

    @property
    def nmm_baset(self) -> ffi.CData:
        return self.__nmm_baset

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_baset_get_abc(self.__nmm_baset)

    def set_lprob(self, nucleotide: bytes, lprob: float) -> None:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        err: int = lib.nmm_baset_set_lprob(self.__nmm_baset, letter, lprob)
        if err != 0:
            nucl = nucleotide.decode()
            raise ValueError(f"Could not set a probability for `{nucl}`.")

    def get_lprob(self, nucleotide: bytes) -> float:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        return lib.nmm_baset_get_lprob(self.__nmm_baset, letter)

    def normalize(self) -> None:
        err: int = lib.nmm_baset_normalize(self.__nmm_baset)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self.__nmm_baset != ffi.NULL:
            lib.nmm_baset_destroy(self.__nmm_baset)


class BaseTable(CBaseTable):
    """
    Nucleotide.

    Parameters
    ----------
    alphabet : Alphabet.
    lprobs : Emission probabilities in log-space.
    """

    def __init__(self, alphabet: Alphabet, lprobs: Dict[bytes, float] = {}):
        super().__init__(alphabet=alphabet)
        self._alphabet = alphabet
        for letter, lprob in lprobs.items():
            self.set_lprob(letter, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet
