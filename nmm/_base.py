from typing import Dict, Optional, Type, TypeVar

from ._alphabet import CAlphabet
from ._ffi import ffi, lib

T = TypeVar("T", bound="BaseTable")


class BaseTable:
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

    def __init__(self, nmm_baset: ffi.CData):
        self._nmm_baset = nmm_baset
        self._alphabet: Optional[CAlphabet] = None
        if self._nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset` is NULL.")

    @classmethod
    def create(cls: Type[T], alphabet: CAlphabet, lprobs: Dict[bytes, float] = {}) -> T:

        nmm_baset = lib.nmm_baset_create(alphabet.imm_abc)
        if nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset_create` failed.")

        baset = cls(nmm_baset)
        baset._alphabet = alphabet
        for letter, lprob in lprobs.items():
            baset.set_lprob(letter, lprob)

        return baset

    @property
    def nmm_baset(self) -> ffi.CData:
        return self._nmm_baset

    @property
    def alphabet(self) -> CAlphabet:
        if self._alphabet is None:
            return CAlphabet(lib.nmm_baset_get_abc(self._nmm_baset))
        return self._alphabet

    def set_lprob(self, nucleotide: bytes, lprob: float) -> None:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        err: int = lib.nmm_baset_set_lprob(self._nmm_baset, letter, lprob)
        if err != 0:
            nucl = nucleotide.decode()
            raise ValueError(f"Could not set a probability for `{nucl}`.")

    def get_lprob(self, nucleotide: bytes) -> float:
        letter = nucleotide
        if len(letter) != 1:
            raise ValueError("Nucleotide must be a single letter.")

        return lib.nmm_baset_get_lprob(self._nmm_baset, letter)

    def normalize(self) -> None:
        err: int = lib.nmm_baset_normalize(self._nmm_baset)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self._nmm_baset != ffi.NULL:
            lib.nmm_baset_destroy(self._nmm_baset)
