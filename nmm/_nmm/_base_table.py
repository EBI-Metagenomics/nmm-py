from typing import Tuple

from .._ffi import ffi, lib
from ._base_alphabet import CBaseAlphabet


class CBaseTable:
    """
    Wrapper around the C implementation of a base table.

    Parameters
    ----------
    nmm_base_table : `<cdata 'struct nmm_base_table *'>`.
        Base table.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_base_table: ffi.CData, base: CBaseAlphabet):
        if nmm_base_table == ffi.NULL:
            raise RuntimeError("`nmm_base_table` is NULL.")
        self._nmm_base_table = nmm_base_table
        self._base = base

    @property
    def alphabet(self) -> CBaseAlphabet:
        return self._base

    @property
    def nmm_base_table(self) -> ffi.CData:
        return self._nmm_base_table

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_base_table_lprob(self._nmm_base_table, nucleotide)

    def __del__(self):
        if self._nmm_base_table != ffi.NULL:
            lib.nmm_base_table_destroy(self._nmm_base_table)


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

    def __init__(self, base: CBaseAlphabet, lprobs: Tuple[float, float, float, float]):
        nmm_base_table = lib.nmm_base_table_create(base.nmm_base_abc, *lprobs)
        super().__init__(nmm_base_table, base)
