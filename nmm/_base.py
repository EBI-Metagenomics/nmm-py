from ._alphabet import Alphabet, CAlphabet
from typing import Dict, Union

from ._ffi import ffi, lib


class CBase:
    """
    Wrapper around the C implementation of a nucleotide.

    Parameters
    ----------
    nmm_base : Union[ffi.CData, None]
        Passing `None` will create a new base at the underlying library level using the `alphabet`
        argument.
    alphabet : Union[CAlphabet, None]
        Passing a `CAlphabet` will create a new base at the underlying library level.
    """

    def __init__(
        self, nmm_base: Union[ffi.CData, None], alphabet: Union[CAlphabet, None]
    ):
        if nmm_base is None:
            if alphabet is None:
                raise ValueError("`alphabet` is `None`")

            self.__cdata = lib.nmm_base_create(alphabet.imm_abc)
        else:
            if alphabet is not None:
                raise ValueError("`alphabet` is not `None`")

            self.__cdata = nmm_base

        if self.__cdata == ffi.NULL:
            raise RuntimeError("`cdata` is NULL.")

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
    """
    Nucleotide.

    Parameters
    ----------
    alphabet : Alphabet.
    lprobs : Emission probabilities in log-space.
    """

    def __init__(self, alphabet: Alphabet, lprobs: Dict[bytes, float] = {}):
        super().__init__(None, alphabet)
        self._alphabet = alphabet
        for letter, lprob in lprobs.items():
            self.set_lprob(letter, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet
