from __future__ import annotations

from typing import Type

from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import Alphabet


class AminoAlphabet(Alphabet):
    """
    Amino acid alphabet is an 20-symbols alphabet.

    Parameters
    ----------
    nmm_amino_abc
        20-symbols alphabet pointer.
    alphabet
        Alphabet.
    """

    def __init__(self, nmm_amino_abc: CData):
        if nmm_amino_abc == ffi.NULL:
            raise RuntimeError("`nmm_amino_abc` is NULL.")
        self._nmm_amino_abc = nmm_amino_abc
        super().__init__(lib.nmm_amino_abc_super(nmm_amino_abc))

    @classmethod
    def create(
        cls: Type[AminoAlphabet], symbols: bytes, any_symbol: bytes
    ) -> AminoAlphabet:
        """
        Create an amino acid alphabet.

        Parameters
        ----------
        symbols
            Set of symbols as an array of bytes.
        any_symbol
            Single-char representing any-symbol.
        """
        if len(any_symbol) != 1:
            raise ValueError("`any_symbol` has length different than 1.")
        return cls(lib.nmm_amino_abc_create(symbols, any_symbol))

    @property
    def nmm_amino_abc(self) -> CData:
        return self._nmm_amino_abc

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class CanonicalAminoAlphabet(AminoAlphabet):
    """
    Canonical amino acid alphabet.

    The canonical symbols are `ACDEFGHIKLMNPQRSTVWY`.
    """

    def __init__(self):
        super().__init__(lib.nmm_amino_abc_create(b"ACDEFGHIKLMNPQRSTVWY", b"X"))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
