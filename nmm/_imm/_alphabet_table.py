from __future__ import annotations

from typing import Sequence, Type

from .._cdata import CData
from .._ffi import ffi, lib
from ._alphabet import Alphabet
from ._lprob import lprob_is_valid


class AlphabetTable:
    """
    Alphabet table of probabilities.

    Parameters
    ----------
    imm_abc_table
        Alphabet table.
    alphabet
        Alphabet.
    """

    def __init__(self, imm_abc_table: CData, alphabet: Alphabet):
        if imm_abc_table == ffi.NULL:
            raise RuntimeError("`imm_abc_table` is NULL.")
        self._imm_abc_table = imm_abc_table
        self._alphabet = alphabet

    @classmethod
    def create(
        cls: Type[AlphabetTable], alphabet: Alphabet, lprobs: Sequence[float]
    ) -> AlphabetTable:
        """
        Clreate an alphabet table of probabilities.

        Parameters
        ----------
        alphabet
            Alphabet.
        lprobs
            Log probability of each nucleotide.
        """
        imm_abc_table = lib.imm_abc_table_create(
            alphabet.imm_abc, ffi.new("double[]", lprobs)
        )
        return cls(imm_abc_table, alphabet)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    @property
    def imm_abc_table(self) -> CData:
        return self._imm_abc_table

    def lprob(self, symbol: bytes) -> float:
        lprob: float = lib.imm_abc_table_lprob(self._imm_abc_table, symbol)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._imm_abc_table != ffi.NULL:
            lib.imm_abc_table_destroy(self._imm_abc_table)
