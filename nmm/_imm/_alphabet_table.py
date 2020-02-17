from typing import Sequence

from .._ffi import ffi, lib
from ._alphabet import CAlphabet
from ._lprob import lprob_is_valid


class CAlphabetTable:
    """
    Wrapper around the C implementation of a alphabet table.

    Parameters
    ----------
    imm_abc_table : `<cdata 'struct imm_abc_table *'>`.
        Alphabet table.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, imm_abc_table: ffi.CData, alphabet: CAlphabet):
        if imm_abc_table == ffi.NULL:
            raise RuntimeError("`imm_abc_table` is NULL.")
        self._imm_abc_table = imm_abc_table
        self._alphabet = alphabet

    @property
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    @property
    def imm_abc_table(self) -> ffi.CData:
        return self._imm_abc_table

    def lprob(self, symbol: bytes) -> float:
        lprob: float = lib.imm_abc_table_lprob(self._imm_abc_table, symbol)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._imm_abc_table != ffi.NULL:
            lib.imm_abc_table_destroy(self._imm_abc_table)


class AlphabetTable(CAlphabetTable):
    """
    Alphabet table of probabilities.

    Parameters
    ----------
    alphabet : `CAlphabet`
        Alphabet.
    lprobs : `Tuple[float, float, float, float]`
        Log probability of each nucleotide.
    """

    def __init__(self, alphabet: CAlphabet, lprobs: Sequence[float]):
        imm_abc_table = lib.imm_abc_table_create(
            alphabet.imm_abc, ffi.new("double[]", lprobs)
        )
        super().__init__(imm_abc_table, alphabet)
