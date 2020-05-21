from __future__ import annotations

from typing import Type

from imm import lprob_is_valid

from ._alphabet import BaseAlphabet
from ._cdata import CData
from ._codon import Codon
from ._ffi import ffi, lib

__all__ = ["CodonProb"]


class CodonProb:
    """
    Codon probabilities.

    Parameters
    ----------
    nmm_codon_lprob
        Codon probabilities pointer.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon_lprob: CData, alphabet: BaseAlphabet):
        if nmm_codon_lprob == ffi.NULL:
            raise RuntimeError("`nmm_codon_lprob` is NULL.")
        self._nmm_codon_lprob = nmm_codon_lprob
        self._alphabet = alphabet

    @classmethod
    def create(cls: Type[CodonProb], alphabet: BaseAlphabet) -> CodonProb:
        """
        Create codon probabilities.

        Parameters
        ----------
        alphabet
            Four-nucleotides alphabet.
        """
        return cls(lib.nmm_codon_lprob_create(alphabet.nmm_base_abc), alphabet)

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def nmm_codon_lprob(self) -> CData:
        return self._nmm_codon_lprob

    def set_lprob(self, codon: Codon, lprob: float):
        if lib.nmm_codon_lprob_set(self._nmm_codon_lprob, codon.nmm_codon, lprob) != 0:
            raise RuntimeError("Could not set codon probability.")

    def get_lprob(self, codon: Codon) -> float:
        lprob: float = lib.nmm_codon_lprob_get(self._nmm_codon_lprob, codon.nmm_codon)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def normalize(self):
        if lib.nmm_codon_lprob_normalize(self._nmm_codon_lprob) != 0:
            raise RuntimeError("Could not normalize.")

    def __del__(self):
        if self._nmm_codon_lprob != ffi.NULL:
            lib.nmm_codon_lprob_destroy(self._nmm_codon_lprob)
