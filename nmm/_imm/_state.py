from typing import Sequence

from .._ffi import ffi, lib
from ._alphabet import Alphabet
from ._lprob import lprob_is_valid
from ._sequence import Sequence as Seq
from ._sequence_table import SequenceTable


class State:
    def __init__(self, imm_state: ffi.CData, alphabet: Alphabet):
        if imm_state == ffi.NULL:
            raise RuntimeError("`imm_state` is NULL.")
        self._imm_state = imm_state
        self._alphabet = alphabet

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    @property
    def imm_state(self) -> ffi.CData:
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

    def lprob(self, sequence: Seq) -> float:
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


class MuteState(State):
    def __init__(self, name: bytes, alphabet: Alphabet):
        """
        Mute state.

        Parameters
        ----------
        name
            State name.
        alphabet
            Alphabet.
        """
        self._imm_mute_state = lib.imm_mute_state_create(name, alphabet.imm_abc)
        if self._imm_mute_state == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")

        super().__init__(lib.imm_state_cast_c(self._imm_mute_state), alphabet)

    def __del__(self):
        if self._imm_mute_state != ffi.NULL:
            lib.imm_mute_state_destroy(self._imm_mute_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class NormalState(State):
    def __init__(self, name: bytes, alphabet: Alphabet, lprobs: Sequence[float]):
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
        state = lib.imm_normal_state_create(name, alphabet.imm_abc, list(lprobs))
        if state == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        self._imm_normal_state = state
        super().__init__(lib.imm_state_cast_c(self._imm_normal_state), alphabet)

    def __del__(self):
        if self._imm_normal_state != ffi.NULL:
            lib.imm_normal_state_destroy(self._imm_normal_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class TableState(State):
    def __init__(self, name: bytes, sequence_table: SequenceTable):
        """
        Parameters
        ----------
        name
            State name.
        sequence_table
            Table of sequence probabilities.
        """
        state = lib.imm_table_state_create(name, sequence_table.imm_seq_table)
        if state == ffi.NULL:
            raise RuntimeError("`imm_table_state_create` failed.")

        self._imm_table_state = state
        alphabet = sequence_table.alphabet
        super().__init__(lib.imm_state_cast_c(self._imm_table_state), alphabet)

    def __del__(self):
        if self._imm_table_state != ffi.NULL:
            lib.imm_table_state_destroy(self._imm_table_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
