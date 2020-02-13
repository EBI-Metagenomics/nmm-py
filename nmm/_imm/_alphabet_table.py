from typing import Sequence

from .._ffi import ffi, lib
from ._alphabet import CAlphabet
from ._lprob import lprob_is_valid


class CAlphabetTable:
    """
    Wrapper around the C implementation of a alphabet table.

    Parameters
    ----------
    imm_abct : `<cdata 'struct imm_abct *'>`.
        Alphabet table.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, imm_abct: ffi.CData, alphabet: CAlphabet):
        if imm_abct == ffi.NULL:
            raise RuntimeError("`imm_abct` is NULL.")
        self._imm_abct = imm_abct
        self._alphabet = alphabet

    @property
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    @property
    def imm_abct(self) -> ffi.CData:
        return self._imm_abct

    def lprob(self, symbol: bytes) -> float:
        lprob: float = lib.imm_abct_lprob(self._imm_abct, symbol)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._imm_abct != ffi.NULL:
            lib.imm_abct_destroy(self._imm_abct)


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
        imm_abct = lib.imm_abct_create(alphabet.imm_abc, ffi.new("double[]", lprobs))
        super().__init__(imm_abct, alphabet)
