from .._ffi import ffi, lib
from .._imm import lprob_is_valid
from ._base_alphabet import CBaseAlphabet
from ._codon import CCodon


class CCodonProb:
    """
    Wrapper around the C implementation of codon probabilities.

    Parameters
    ----------
    nmm_codon_lprob : `<cdata 'struct nmm_codon_lprob *'>`.
        Codon probabilities pointer.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon_lprob: ffi.CData, base: CBaseAlphabet):
        if nmm_codon_lprob == ffi.NULL:
            raise RuntimeError("`nmm_codon_lprob` is NULL.")
        self._nmm_codon_lprob = nmm_codon_lprob
        self._base = base

    @property
    def base(self) -> CBaseAlphabet:
        return self._base

    @property
    def nmm_codon_lprob(self) -> ffi.CData:
        return self._nmm_codon_lprob

    def set_lprob(self, codon: CCodon, lprob: float):
        if lib.nmm_codon_lprob_set(self._nmm_codon_lprob, codon.nmm_codon, lprob) != 0:
            raise RuntimeError("Could not set codon probability.")

    def get_lprob(self, codon: CCodon) -> float:
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


class CodonProb(CCodonProb):
    """
    Codon probabilities.

    Parameters
    ----------
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, base: CBaseAlphabet):
        super().__init__(lib.nmm_codon_lprob_create(base.nmm_base_abc), base)
