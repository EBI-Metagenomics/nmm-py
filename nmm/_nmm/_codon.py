from .._ffi import ffi, lib
from ._base import CBase


class CCodon:
    """
    Wrapper around the C implementation of codon.

    Codon is a triplet of symbols from a given alphabet.

    Parameters
    ----------
    nmm_codon : `<cdata 'struct nmm_codon *'>`.
        Codon pointer.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon: ffi.CData, base: CBase):
        if nmm_codon == ffi.NULL:
            raise RuntimeError("`nmm_codon` is NULL.")
        self._nmm_codon = nmm_codon
        self._base = base

    @property
    def base(self) -> CBase:
        return self._base

    @property
    def symbols(self) -> bytes:
        triplet = lib.nmm_codon_get_triplet(self._nmm_codon)
        return triplet.a + triplet.b + triplet.c

    @symbols.setter
    def symbols(self, symbols: bytes):
        if len(symbols) != 3:
            raise ValueError("Symbols length must be three.")

        triplet = {"a": symbols[0:1], "b": symbols[1:2], "c": symbols[2:3]}
        if lib.nmm_codon_set_triplet(self._nmm_codon, triplet) != 0:
            raise ValueError("Could not set codon.")

    @property
    def nmm_codon(self) -> ffi.CData:
        return self._nmm_codon

    def __del__(self):
        if self._nmm_codon != ffi.NULL:
            lib.nmm_codon_destroy(self._nmm_codon)

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Codon(CCodon):
    """
    Codon is a sequence of three symbols from a four-nucleotides alphabet.

    Parameters
    ----------
    symbols : bytes
        Sequence of four symbols.
    base : `CBase`
        Four-nucleotides alphabet.
    """

    def __init__(self, symbols: bytes, base: CBase):
        super().__init__(lib.nmm_codon_create(base.nmm_base), base)
        self.symbols = symbols

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
