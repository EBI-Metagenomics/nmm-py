from ._alphabet import CAlphabet
from ._ffi import ffi, lib


class CBase:
    """
    Wrapper around the C implementation of a base (four-nucleotides alphabet).

    Parameters
    ----------
    nmm_base : `<cdata 'struct nmm_base *'>`.
        Four-nucleotides alphabet pointer.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, nmm_base: ffi.CData, alphabet: CAlphabet):
        if nmm_base == ffi.NULL:
            raise RuntimeError("`nmm_base` is NULL.")
        if lib.nmm_base_get_abc(nmm_base) != alphabet.imm_abc:
            raise ValueError("Alphabets must be the same.")
        self._nmm_base = nmm_base
        self._alphabet = alphabet

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def nmm_base(self) -> ffi.CData:
        return self._nmm_base

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_abc_symbols(lib.nmm_base_get_abc(self._nmm_base)))

    def __del__(self):
        if self._nmm_base != ffi.NULL:
            lib.nmm_base_destroy(self._nmm_base)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Base(CBase):
    """
    Base is a four-nucleotides alphabet.

    Parameters
    ----------
    alphabet : `Alphabet`
        Four-nucleotides alphabet.
    """

    def __init__(self, alphabet: CAlphabet):
        super().__init__(lib.nmm_base_create(alphabet.imm_abc), alphabet)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
