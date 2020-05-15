from __future__ import annotations

from typing import Type

from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import Alphabet


class BaseAlphabet(Alphabet):
    """
    Base alphabet is a four-nucleotides alphabet.

    Parameters
    ----------
    nmm_base_abc
        Four-nucleotides alphabet pointer.
    alphabet
        Alphabet.
    """

    def __init__(self, nmm_base_abc: CData):
        if nmm_base_abc == ffi.NULL:
            raise RuntimeError("`nmm_base_abc` is NULL.")
        self._nmm_base_abc = nmm_base_abc
        super().__init__(lib.nmm_base_abc_super(nmm_base_abc))

    @classmethod
    def create(
        cls: Type[BaseAlphabet], symbols: bytes, any_symbol: bytes
    ) -> BaseAlphabet:
        """
        Create a base alphabet.

        Parameters
        ----------
        symbols
            Set of symbols as an array of bytes.
        any_symbol
            Single-char representing any-symbol.
        """
        if len(any_symbol) != 1:
            raise ValueError("`any_symbol` has length different than 1.")
        return cls(lib.nmm_base_abc_create(symbols, any_symbol))

    @property
    def nmm_base_abc(self) -> CData:
        return self._nmm_base_abc

    def __del__(self):
        if self._nmm_base_abc != ffi.NULL:
            lib.nmm_base_abc_destroy(self._nmm_base_abc)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class RNAAlphabet(BaseAlphabet):
    """
    RNA alphabet.
    """

    def __init__(self):
        abc = Alphabet.create(b"ACGU", b"X")
        super().__init__(lib.nmm_base_abc_create(abc.imm_abc), abc)


class DNAAlphabet(BaseAlphabet):
    """
    DNA alphabet.
    """

    def __init__(self):
        abc = Alphabet.create(b"ACGT", b"X")
        super().__init__(lib.nmm_base_abc_create(abc.imm_abc), abc)
