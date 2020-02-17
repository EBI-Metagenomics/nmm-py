from .._ffi import ffi, lib
from .._imm import CAlphabet


class CBase:
    """
    Wrapper around the C implementation of a base (four-nucleotides alphabet).

    Parameters
    ----------
    nmm_base_abc : `<cdata 'struct nmm_base_abc *'>`.
        Four-nucleotides alphabet pointer.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, nmm_base_abc: ffi.CData, alphabet: CAlphabet):
        if nmm_base_abc == ffi.NULL:
            raise RuntimeError("`nmm_base_abc` is NULL.")
        if lib.nmm_base_abc_cast(nmm_base_abc) != alphabet.imm_abc:
            raise ValueError("Alphabets must be the same.")
        self._nmm_base_abc = nmm_base_abc
        self._alphabet = alphabet

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def nmm_base_abc(self) -> ffi.CData:
        return self._nmm_base_abc

    @property
    def symbols(self) -> bytes:
        return ffi.string(
            lib.imm_abc_symbols(lib.nmm_base_abc_cast(self._nmm_base_abc))
        )

    def __del__(self):
        if self._nmm_base_abc != ffi.NULL:
            lib.nmm_base_abc_destroy(self._nmm_base_abc)

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
        super().__init__(lib.nmm_base_abc_create(alphabet.imm_abc), alphabet)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
