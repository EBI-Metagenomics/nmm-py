from __future__ import annotations
from typing import Type
from .._ffi import ffi, lib


class Alphabet:
    """
    Set of symbols for Markov Models.

    Parameters
    ----------
    imm_abc
        Pointer.
    """

    def __init__(self, imm_abc: ffi.CData):
        super().__init__()
        if imm_abc == ffi.NULL:
            raise RuntimeError("`imm_abc` is NULL.")
        self._imm_abc = imm_abc

    @classmethod
    def create(cls: Type[Alphabet], symbols: bytes, any_symbol: bytes) -> Alphabet:
        """
        Create an alphabet.

        Parameters
        ----------
        symbols
            Set of symbols as an array of bytes.
        any_symbol
            Single-char representing any-symbol.
        """
        if len(any_symbol) != 1:
            raise ValueError("`any_symbol` has length different than 1.")
        super().__init__(lib.imm_abc_create(symbols, any_symbol))
        return cls(lib.imm_abc_create(symbols, any_symbol))

    @property
    def imm_abc(self) -> ffi.CData:
        return self._imm_abc

    @property
    def length(self) -> int:
        return lib.imm_abc_length(self._imm_abc)

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_abc_symbols(self._imm_abc))

    def has_symbol(self, symbol_id: bytes) -> bool:
        return lib.imm_abc_has_symbol(self._imm_abc, symbol_id)

    def symbol_idx(self, symbol_id: bytes) -> int:
        return lib.imm_abc_symbol_idx(self._imm_abc, symbol_id)

    def symbol_id(self, symbol_idx: int) -> bytes:
        return lib.imm_abc_symbol_id(self._imm_abc, symbol_idx)

    def __del__(self):
        if self._imm_abc != ffi.NULL:
            lib.imm_abc_destroy(self._imm_abc)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
