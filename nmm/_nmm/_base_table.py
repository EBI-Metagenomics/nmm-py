from typing import Tuple

from .._ffi import ffi, lib
from ._base import CBase


class CBaseTable:
    """
    Wrapper around the C implementation of a base table.

    Parameters
    ----------
    nmm_baset : `<cdata 'struct nmm_baset *'>`.
        Base table.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_baset: ffi.CData, base: CBase):
        if nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset` is NULL.")
        self._nmm_baset = nmm_baset
        self._base = base

    @property
    def base(self) -> CBase:
        return self._base

    @property
    def nmm_baset(self) -> ffi.CData:
        return self._nmm_baset

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_baset_lprob(self._nmm_baset, nucleotide)

    def __del__(self):
        if self._nmm_baset != ffi.NULL:
            lib.nmm_baset_destroy(self._nmm_baset)


class BaseTable(CBaseTable):
    """
    Base table of probabilities.

    Parameters
    ----------
    base : `CBase`
        Four-nucleotides alphabet.
    lprobs : `Tuple[float, float, float, float]`
        Log probability of each nucleotide.
    """

    def __init__(self, base: CBase, lprobs: Tuple[float, float, float, float]):
        nmm_baset = lib.nmm_baset_create(base.nmm_base, *lprobs)
        super().__init__(nmm_baset, base)
