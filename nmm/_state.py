from typing import Sequence

from ._alphabet import CAlphabet
from ._base_table import BaseTable
from ._codon_table import CodonTable
from ._codon import Codon
from ._ffi import ffi, lib
from ._lprob import lprob_is_valid
from ._sequence import CSequence
from ._sequence_table import SequenceTable


class CState:
    def __init__(self, imm_state: ffi.CData):
        if imm_state == ffi.NULL:
            raise RuntimeError("`imm_state` is NULL.")
        self._imm_state = imm_state

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

    def lprob(self, seq: CSequence) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        seq : `CSequence`
            Sequence.
        """
        lprob: float = lib.imm_state_lprob(self._imm_state, seq.imm_seq)
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


class MuteState(CState):
    def __init__(self, name: bytes, alphabet: CAlphabet):
        """
        Mute state.

        Parameters
        ----------
        name : bytes
            State name.
        alphabet : `CAlphabet`
            Alphabet.
        """
        self._imm_mute_state = lib.imm_mute_state_create(name, alphabet.imm_abc)
        if self._imm_mute_state == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")

        super().__init__(lib.imm_state_cast_c(self._imm_mute_state))

    def __del__(self):
        if self._imm_mute_state != ffi.NULL:
            lib.imm_mute_state_destroy(self._imm_mute_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class NormalState(CState):
    def __init__(self, name: bytes, alphabet: CAlphabet, lprobs: Sequence[float]):
        """
        Normal state.

        Parameters
        ----------
        name : bytes
            State name.
        alphabet : `CAlphabet`
            Alphabet.
        lprobs : `typing.Sequence[float]`
            Emission probabilities in log-space for each alphabet letter.
        """
        state = lib.imm_normal_state_create(name, alphabet.imm_abc, list(lprobs))
        if state == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        self._imm_normal_state = state
        super().__init__(lib.imm_state_cast_c(self._imm_normal_state))

    def __del__(self):
        if self._imm_normal_state != ffi.NULL:
            lib.imm_normal_state_destroy(self._imm_normal_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class TableState(CState):
    def __init__(self, name: bytes, sequence_table: SequenceTable):
        """
        Parameters
        ----------
        name : bytes
            State name.
        sequence_table : `SequenceTable`
            Table of sequence probabilities.
        """
        state = lib.imm_table_state_create(name, sequence_table.imm_seq_table)
        if state == ffi.NULL:
            raise RuntimeError("`imm_table_state_create` failed.")

        self._imm_table_state = state
        super().__init__(lib.imm_state_cast_c(self._imm_table_state))

    def __del__(self):
        if self._imm_table_state != ffi.NULL:
            lib.imm_table_state_destroy(self._imm_table_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class FrameState(CState):
    def __init__(
        self, name: bytes, baset: BaseTable, codont: CodonTable, epsilon: float
    ):
        """
        Parameters
        ----------
        name : bytes
            State name.
        baset : `BaseTable`
            Base table of probabilities.
        codont : `CodonTable`
            Codon table of probabilities.
        epsilon : float
            Epsilon.
        """
        state = lib.nmm_frame_state_create(
            name, baset.nmm_baset, codont.nmm_codont, epsilon
        )
        if state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state_create` failed.")

        self._baset = baset
        self._codont = codont
        self._epsilon = epsilon
        self._nmm_frame_state = state
        super().__init__(lib.imm_state_cast_c(self._nmm_frame_state))

    def decode(self, seq: CSequence, codon: Codon) -> float:
        state = self._nmm_frame_state
        return lib.nmm_frame_state_decode(state, seq.imm_seq, codon.nmm_codon)

    def __del__(self):
        if self._nmm_frame_state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._nmm_frame_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
