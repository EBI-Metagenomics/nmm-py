from __future__ import annotations
from enum import Enum
from typing import Generic, Iterable, TypeVar, Type

from .._cdata import CData
from .._ffi import ffi, lib
from ._alphabet import Alphabet
from ._lprob import lprob_is_valid
from ._sequence import Sequence
from ._sequence_table import SequenceTable

T = TypeVar("T", bound=Alphabet)


class StateType(Enum):
    MUTE = 0x00
    NORMAL = 0x01
    TABLE = 0x02


class State(Generic[T]):
    def __init__(self, imm_state: CData, alphabet: T):
        if imm_state == ffi.NULL:
            raise RuntimeError("`imm_state` is NULL.")
        self._imm_state = imm_state
        self._alphabet = alphabet

    @property
    def alphabet(self) -> T:
        return self._alphabet

    @property
    def imm_state(self) -> CData:
        return self._imm_state

    @property
    def name(self) -> bytes:
        return ffi.string(lib.imm_state_get_name(self._imm_state))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self._imm_state)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self._imm_state)

    def lprob(self, sequence: Sequence) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        sequence
            Sequence.
        """
        lprob: float = lib.imm_state_lprob(self._imm_state, sequence.imm_seq)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __str__(self) -> str:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self._imm_state == ffi.NULL:
            raise RuntimeError("State has failed to initialize.")
        return f"{self.name.decode()}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class MuteState(State[T]):
    def __init__(self, imm_mute_state: CData, alphabet: T):
        """
        Mute state.

        Parameters
        ----------
        imm_mute_state
            State pointer.
        alphabet
            Alphabet.
        """
        if imm_mute_state == ffi.NULL:
            raise RuntimeError("`imm_mute_state` is NULL.")
        self._imm_mute_state = imm_mute_state
        super().__init__(lib.imm_mute_state_super(self._imm_mute_state), alphabet)

    @classmethod
    def create(cls: Type[MuteState], name: bytes, alphabet: T) -> MuteState:
        """
        Mute state.

        Parameters
        ----------
        name
            State name.
        alphabet
            Alphabet.
        """
        imm_mute_state = lib.imm_mute_state_create(name, alphabet.imm_abc)
        return MuteState(imm_mute_state, alphabet)

    def __del__(self):
        if self._imm_mute_state != ffi.NULL:
            lib.imm_mute_state_destroy(self._imm_mute_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class NormalState(State[T]):
    def __init__(self, imm_normal_state: CData, alphabet: T):
        """
        Normal state.

        Parameters
        ----------
        imm_normal_state
            State pointer.
        alphabet
            Alphabet.
        """
        if imm_normal_state == ffi.NULL:
            raise RuntimeError("`imm_normal_state` is NULL.")
        self._imm_normal_state = imm_normal_state
        super().__init__(lib.imm_normal_state_super(self._imm_normal_state), alphabet)

    @classmethod
    def create(
        cls: Type[NormalState], name: bytes, alphabet: T, lprobs: Iterable[float]
    ) -> NormalState:
        """
        Normal state.

        Parameters
        ----------
        name
            State name.
        alphabet
            Alphabet.
        lprobs
            Emission probabilities in log-space for each alphabet letter.
        """
        ptr = lib.imm_normal_state_create(name, alphabet.imm_abc, list(lprobs))
        return NormalState(ptr, alphabet)

    def __del__(self):
        if self._imm_normal_state != ffi.NULL:
            lib.imm_normal_state_destroy(self._imm_normal_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class TableState(State[T]):
    def __init__(self, imm_table_state: CData, alphabet: T):
        """
        Table state.

        Parameters
        ----------
        imm_table_state
            State pointer.
        alphabet
            Alphabet.
        """
        if imm_table_state == ffi.NULL:
            raise RuntimeError("`imm_table_state` is NULL.")
        self._imm_table_state = imm_table_state
        super().__init__(lib.imm_table_state_super(self._imm_table_state), alphabet)

    @classmethod
    def create(
        cls: Type[TableState], name: bytes, sequence_table: SequenceTable
    ) -> TableState:
        """
        Table state.

        Parameters
        ----------
        name
            State name.
        sequence_table
            Table of sequence probabilities.
        """
        ptr = lib.imm_table_state_create(name, sequence_table.imm_seq_table)
        return TableState(ptr, sequence_table.alphabet)

    def __del__(self):
        if self._imm_table_state != ffi.NULL:
            lib.imm_table_state_destroy(self._imm_table_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


def wrap_imm_state(imm_state: CData) -> State:
    state_type: int = lib.imm_state_type_id(imm_state)
    if state_type == StateType.MUTE:
        pass
    elif state_type == StateType.NORMAL:
        pass
    elif state_type == StateType.TABLE:
        pass
    raise ValueError(f"Unknown state type: {state_type}.")
