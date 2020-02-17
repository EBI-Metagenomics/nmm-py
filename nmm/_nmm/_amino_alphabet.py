from .._ffi import ffi, lib
from .._imm import CAlphabet


class CAminoAlphabet:
    """
    Wrapper around the C implementation of a amino alphabet.

    Parameters
    ----------
    nmm_amino_abc : `<cdata 'struct nmm_amino_abc *'>`.
        20 symbols alphabet pointer.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, nmm_amino_abc: ffi.CData, alphabet: CAlphabet):
        if nmm_amino_abc == ffi.NULL:
            raise RuntimeError("`nmm_amino_abc` is NULL.")
        if lib.nmm_amino_abc_cast(nmm_amino_abc) != alphabet.imm_abc:
            raise ValueError("Alphabets must be the same.")
        self._nmm_amino_abc = nmm_amino_abc
        self._alphabet = alphabet

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def nmm_amino_abc(self) -> ffi.CData:
        return self._nmm_amino_abc

    @property
    def symbols(self) -> bytes:
        return ffi.string(
            lib.imm_abc_symbols(lib.nmm_amino_abc_cast(self._nmm_amino_abc))
        )

    def __del__(self):
        if self._nmm_amino_abc != ffi.NULL:
            lib.nmm_amino_abc_destroy(self._nmm_amino_abc)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class AminoAlphabet(CAminoAlphabet):
    """
    Amino alphabet is an 20-symbols alphabet.

    Parameters
    ----------
    alphabet : `Alphabet`
        20 symbols alphabet.
    """

    def __init__(self, alphabet: CAlphabet):
        super().__init__(lib.nmm_amino_abc_create(alphabet.imm_abc), alphabet)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"