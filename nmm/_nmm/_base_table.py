from __future__ import annotations

from typing import Tuple, Type

from .._cdata import CData
from .._ffi import ffi, lib
from ._base_alphabet import BaseAlphabet


class BaseTable:
    """
    Base table of probabilities.

    Parameters
    ----------
    nmm_base_table
        Base table.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_base_table: CData, alphabet: BaseAlphabet):
        if nmm_base_table == ffi.NULL:
            raise RuntimeError("`nmm_base_table` is NULL.")
        self._nmm_base_table = nmm_base_table
        self._alphabet = alphabet

    @classmethod
    def create(
        cls: Type[BaseTable],
        alphabet: BaseAlphabet,
        lprobs: Tuple[float, float, float, float],
    ) -> BaseTable:
        """
        Create base table of probabilities.

        Parameters
        ----------
        alphabet
            Four-nucleotides alphabet.
        lprobs
            Log probability of each nucleotide.
        """
        nmm_base_table = lib.nmm_base_table_create(alphabet.nmm_base_abc, *lprobs)
        return cls(nmm_base_table, alphabet)

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def nmm_base_table(self) -> CData:
        return self._nmm_base_table

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_base_table_lprob(self._nmm_base_table, nucleotide)

    def __del__(self):
        if self._nmm_base_table != ffi.NULL:
            lib.nmm_base_table_destroy(self._nmm_base_table)
