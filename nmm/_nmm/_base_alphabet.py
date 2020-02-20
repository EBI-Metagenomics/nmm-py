from __future__ import annotations

from typing import Type

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

    def __init__(self, nmm_base_abc: ffi.CData, alphabet: Alphabet):
        if nmm_base_abc == ffi.NULL:
            raise RuntimeError("`nmm_base_abc` is NULL.")
        if lib.nmm_base_abc_cast(nmm_base_abc) != alphabet.imm_abc:
            raise ValueError("Alphabets must be the same.")
        self._nmm_base_abc = nmm_base_abc
        self._alphabet = alphabet
        super().__init__(lib.nmm_base_abc_cast(nmm_base_abc))

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
        abc = Alphabet.create(symbols, any_symbol)
        return cls(lib.nmm_base_abc_create(abc.imm_abc), abc)

    @property
    def nmm_base_abc(self) -> ffi.CData:
        return self._nmm_base_abc

    def __del__(self):
        if self._nmm_base_abc != ffi.NULL:
            lib.nmm_base_abc_destroy(self._nmm_base_abc)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
