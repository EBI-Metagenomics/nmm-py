from ._base import CBase
from ._codon import CCodon
from ._ffi import ffi, lib
from ._lprob import lprob_is_valid


class CCodonProb:
    """
    Wrapper around the C implementation of a codon probability.

    Parameters
    ----------
    nmm_codonp : `<cdata 'struct nmm_codonp *'>`.
    """

    def __init__(self, nmm_codonp: ffi.CData):
        super().__init__()
        if nmm_codonp == ffi.NULL:
            raise RuntimeError("`nmm_codonp` is NULL.")
        self._nmm_codonp = nmm_codonp

    @property
    def nmm_codonp(self) -> ffi.CData:
        return self._nmm_codonp

    def set_lprob(self, codon: CCodon, lprob: float):
        if lib.nmm_codonp_set_lprob(self._nmm_codonp, codon.nmm_codon, lprob) != 0:
            raise RuntimeError("Could not set codon probability.")

    def get_lprob(self, codon: CCodon) -> float:
        lprob: float = lib.nmm_codonp_get_lprob(self._nmm_codonp, codon.nmm_codon)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def normalize(self):
        if lib.nmm_codonp_normalize(self._nmm_codonp) != 0:
            raise RuntimeError("Could not normalize.")

    def __del__(self):
        if self._nmm_codonp != ffi.NULL:
            lib.nmm_codonp_destroy(self._nmm_codonp)


class CodonProb(CCodonProb):
    """
    Codon probabilities.

    Parameters
    ----------
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, base: CBase):
        self._calphabet = base
        nmm_codonp = lib.nmm_codonp_create(base.nmm_base)
        super().__init__(nmm_codonp)
