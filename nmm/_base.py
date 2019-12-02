from typing import Dict, Optional, Type, TypeVar, Union

from ._alphabet import CAlphabet
from ._ffi import ffi, lib


class Base:
    def __init__(self, base: Union[bytes, str]):
        if isinstance(base, str):
            base = base.encode()

        if len(base) != 1:
            raise ValueError("Base must be a single letter.")

        self._base = base

    def __bytes__(self) -> bytes:
        return self._base

    def __str__(self) -> str:
        return self._base.decode()

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{self._base.decode()}>"


T = TypeVar("T", bound="BaseTable")


class BaseTable:
    """
    Wrapper around the C implementation of a base table.

    Parameters
    ----------
    nmm_baset : CData
        Base table.
    """

    def __init__(self, nmm_baset: ffi.CData):
        self._nmm_baset = nmm_baset
        self._alphabet: Optional[CAlphabet] = None
        if self._nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset` is NULL.")

    @classmethod
    def create(cls: Type[T], alphabet: CAlphabet, lprobs: Dict[Base, float] = {}) -> T:

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

    def set_lprob(self, base: Base, lprob: float) -> None:
        err: int = lib.nmm_baset_set_lprob(self._nmm_baset, bytes(base), lprob)
        if err != 0:
            nucl = str(base)
            raise ValueError(f"Could not set a probability for `{nucl}`.")

    def get_lprob(self, base: Base) -> float:
        return lib.nmm_baset_get_lprob(self._nmm_baset, bytes(base))

    def normalize(self) -> None:
        err: int = lib.nmm_baset_normalize(self._nmm_baset)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self._nmm_baset != ffi.NULL:
            lib.nmm_baset_destroy(self._nmm_baset)
