from ._codon import CCodon
from ._codon_prob import CCodonProb
from ._ffi import ffi, lib
from ._lprob import lprob_is_valid


class CCodonTable:
    """
    Wrapper around the C implementation of a codon table.

    Parameters
    ----------
    nmm_codont : `<cdata 'struct nmm_codont *'>`.
    """

    def __init__(self, nmm_codont: ffi.CData):
        if nmm_codont == ffi.NULL:
            raise RuntimeError("`nmm_codont` is NULL.")
        self._nmm_codont = nmm_codont

    @property
    def nmm_codont(self) -> ffi.CData:
        return self._nmm_codont

    def lprob(self, codon: CCodon) -> float:
        lprob: float = lib.nmm_codont_lprob(self._nmm_codont, codon.nmm_codon)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._nmm_codont != ffi.NULL:
            lib.nmm_codont_destroy(self._nmm_codont)


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
        self._codonp = codonp
        super().__init__(lib.nmm_codont_create(codonp.nmm_codonp))
