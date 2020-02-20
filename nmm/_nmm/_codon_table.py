from __future__ import annotations

from typing import Type

from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import lprob_is_valid
from ._base_alphabet import BaseAlphabet
from ._codon import Codon
from ._codon_prob import CodonProb


class CodonTable:
    """
    Codon table.

    Compute marginal and non-marginal codon probabilities.

    Parameters
    ----------
    nmm_codon_table
        Codon table.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon_table: CData, alphabet: BaseAlphabet):
        if nmm_codon_table == ffi.NULL:
            raise RuntimeError("`nmm_codon_table` is NULL.")
        self._nmm_codon_table = nmm_codon_table
        self._alphabet = alphabet

    @classmethod
    def create(
        cls: Type[CodonTable],
        codonp: CodonProb,
    ) -> CodonTable:
        """
        Create a codon table.

        Parameters
        ----------
        codonp
            Non-marginal codon probabilities.
        """
        return cls(lib.nmm_codon_table_create(codonp.nmm_codon_lprob), codonp.alphabet)

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def nmm_codon_table(self) -> CData:
        return self._nmm_codon_table

    def lprob(self, codon: Codon) -> float:
        lprob: float = lib.nmm_codon_table_lprob(self._nmm_codon_table, codon.nmm_codon)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._nmm_codon_table != ffi.NULL:
            lib.nmm_codon_table_destroy(self._nmm_codon_table)
