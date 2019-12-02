from abc import ABC, abstractmethod
from typing import Type, TypeVar

from ._ffi import ffi, lib


class AlphabetBase(ABC):
    @property
    @abstractmethod
    def length(self) -> int:
        raise NotImplementedError()

    @property
    @abstractmethod
    def symbols(self) -> bytes:
        raise NotImplementedError()

    @abstractmethod
    def has_symbol(self, symbol_id: bytes) -> bool:
        del symbol_id
        raise NotImplementedError()

    @abstractmethod
    def symbol_idx(self, symbol_id: bytes) -> int:
        del symbol_id
        raise NotImplementedError()

    @abstractmethod
    def symbol_id(self, symbol_idx: int) -> bytes:
        del symbol_idx
        raise NotImplementedError()

    def __str__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{symbols}}}"

    def __repr__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{self.__class__.__name__}:{symbols}}}"


T = TypeVar("T", bound="CAlphabet")


class CAlphabet(AlphabetBase):
    """
    Wrapper around the C implementation of alphabet set.

    Parameters
    ----------
    imm_abc : `<cdata 'struct imm_abc *'>`.
    """

    def __init__(self, imm_abc: ffi.CData):
        self._imm_abc = imm_abc

    @classmethod
    def clone_from_imm_abc(cls: Type[T], imm_abc: ffi.CData) -> T:
        t = lib.imm_abc_clone(imm_abc)
        if t == ffi.NULL:
            raise RuntimeError("`imm_abc_clone` failed.")
        return cls(t)

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
        return lib.imm_abc_has_symbol(self._imm_abc, symbol_id) == 1

    def symbol_idx(self, symbol_id: bytes) -> int:
        return lib.imm_abc_symbol_idx(self._imm_abc, symbol_id)

    def symbol_id(self, symbol_idx: int) -> bytes:
        return lib.imm_abc_symbol_id(self._imm_abc, symbol_idx)

    def __del__(self):
        if self._imm_abc != ffi.NULL:
            lib.imm_abc_destroy(self._imm_abc)

    def __repr__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{self.__class__.__name__}:{symbols}}}"


class Alphabet(CAlphabet):
    """
    Alphabet set for Markov models.

    Parameters
    ----------
    symbols : bytes
        Set of symbols as an array of bytes.
    """

    def __init__(self, symbols: bytes, any_symbol: bytes = b"X"):
        if len(any_symbol) != 1:
            raise ValueError("`any_symbol` has length different than 1.")

        imm_abc = lib.imm_abc_create(symbols, any_symbol)
        if imm_abc == ffi.NULL:
            raise RuntimeError("`imm_abc_create` failed.")

        super().__init__(imm_abc)

    def __repr__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{self.__class__.__name__}:{symbols}}}"
