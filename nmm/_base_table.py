from typing import Tuple
from ._base import CBase
from ._ffi import ffi, lib


class CBaseTable:
    """
    Wrapper around the C implementation of a base table.

    Parameters
    ----------
    nmm_baset : `<cdata 'struct nmm_baset *'>`.
    """

    def __init__(self, nmm_baset: ffi.CData):
        if nmm_baset == ffi.NULL:
            raise RuntimeError("`nmm_baset` is NULL.")
        self._nmm_baset = nmm_baset

    @property
    def nmm_baset(self) -> ffi.CData:
        return self._nmm_baset

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_baset_lprob(self._nmm_baset, nucleotide)

    def __del__(self):
        if self._nmm_baset != ffi.NULL:
            lib.nmm_baset_destroy(self._nmm_baset)


class BaseTable(CBaseTable):
    def __init__(self, base: CBase, lprobs: Tuple[float, float, float, float]):
        self._base = base
        nmm_baset = lib.nmm_baset_create(base.nmm_base, *lprobs)
        super().__init__(nmm_baset)
