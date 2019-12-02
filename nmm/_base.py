from typing import Dict, Type, TypeVar, Union, Optional

from ._alphabet import CAlphabet
from ._ffi import ffi, lib


class Base:
    """
    Base is a nucleotide letter.

    base : Union[bytes, str, int]
        A single letter.
    """

    def __init__(self, base: Union[bytes, str, int]):
        if isinstance(base, str):
            base = base.encode()

        if isinstance(base, int):
            base = bytes([base])

        if len(base) != 1:
            raise ValueError("Base must be a single letter.")

        self._base = base

    def __eq__(self, another):
        return bytes(self) == bytes(another)

    def __hash__(self):
        return hash(bytes(self))

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
        if nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset` is NULL.")
        self._nmm_baset = nmm_baset
        self._alphabet: Optional[CAlphabet] = None

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
            imm_abc = lib.nmm_baset_get_abc(self._nmm_baset)
            return CAlphabet.clone_from_imm_abc(imm_abc)
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
