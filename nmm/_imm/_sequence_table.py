from __future__ import annotations

from typing import Type, TypeVar, Generic

from .._cdata import CData
from .._ffi import ffi, lib
from ._alphabet import Alphabet
from ._lprob import lprob_is_valid
from ._sequence import Sequence

T = TypeVar("T", bound=Alphabet)


class SequenceTable(Generic[T]):
    """
    Table of sequence probabilities.

    Parameters
    ----------
    imm_seq_table
        Pointer.
    alphabet
        Alphabet.
    """

    def __init__(self, imm_seq_table: CData, alphabet: T):
        if imm_seq_table == ffi.NULL:
            raise RuntimeError("`imm_seq_table` is NULL.")
        self._imm_seq_table = imm_seq_table
        self._alphabet = alphabet

    @classmethod
    def create(cls: Type[SequenceTable], alphabet: T) -> SequenceTable:
        """
        Create a table of sequence probabilities.

        Parameters
        ----------
        alphabet
            Alphabet.
        """
        return cls(lib.imm_seq_table_create(alphabet.imm_abc), alphabet)

    @property
    def alphabet(self) -> T:
        return self._alphabet

    @property
    def imm_seq_table(self) -> CData:
        return self._imm_seq_table

    @property
    def length(self) -> int:
        imm_abc = lib.imm_seq_table_get_abc(self._imm_seq_table)
        return lib.imm_abc_length(imm_abc)

    @property
    def symbols(self) -> bytes:
        imm_abc = lib.imm_seq_table_get_abc(self._imm_seq_table)
        return ffi.string(lib.imm_abc_string(imm_abc))

    def add(self, sequence: Sequence, lprob: float):
        if lib.imm_seq_table_add(self._imm_seq_table, sequence.imm_seq, lprob) != 0:
            raise RuntimeError("Could not add sequence.")

    def normalize(self):
        if lib.imm_seq_table_normalize(self._imm_seq_table) != 0:
            raise RuntimeError("Could not normalize it.")

    def lprob(self, sequence: Sequence) -> float:
        lprob: float = lib.imm_seq_table_lprob(self._imm_seq_table, sequence.imm_seq)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._imm_seq_table != ffi.NULL:
            lib.imm_seq_table_destroy(self._imm_seq_table)

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
