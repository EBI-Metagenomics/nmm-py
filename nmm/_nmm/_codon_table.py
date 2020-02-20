from .._ffi import ffi, lib
from .._imm import lprob_is_valid
from ._base_alphabet import BaseAlphabet
from ._codon import Codon
from .._cdata import CData
from ._codon_prob import CCodonProb


class CCodonTable:
    """
    Wrapper around the C implementation of a codon table.

    Parameters
    ----------
    nmm_codon_table : `<cdata 'struct nmm_codon_table *'>`.
        Codon table.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon_table: CData, base: BaseAlphabet):
        if nmm_codon_table == ffi.NULL:
            raise RuntimeError("`nmm_codon_table` is NULL.")
        self._nmm_codon_table = nmm_codon_table
        self._base = base

    @property
    def base(self) -> BaseAlphabet:
        return self._base

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


class CodonTable(CCodonTable):
    """
    Codon table.

    Compute marginal and non-marginal codon probabilities.

    Parameters
    ----------
    codonp : `CCodonProb`
        Non-marginal codon probabilities.
    """

    def __init__(self, codonp: CCodonProb):
        super().__init__(
            lib.nmm_codon_table_create(codonp.nmm_codon_lprob), codonp.base
        )
