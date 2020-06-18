from __future__ import annotations

from enum import Enum
from typing import Tuple, Type, TypeVar

from imm import Alphabet, Sequence, State

from ._alphabet import BaseAlphabet
from ._cdata import CData
from ._codon import Codon
from ._codon_prob import CodonProb
from ._ffi import ffi, lib
from ._table import BaseTable, CodonTable

__all__ = [
    "CodonState",
    "FrameState",
    "StateType",
]

T = TypeVar("T", bound=Alphabet)


class StateType(Enum):
    CODON = 0x10
    FRAME = 0x11


class FrameState(State[BaseAlphabet]):
    def __init__(
        self, nmm_frame_state: CData, baset: BaseTable, codont: CodonTable,
    ):
        """
        Frame state.

        Parameters
        ----------
        nmm_frame_state
            State pointer.
        baset
            Base table of probabilities.
        codont
            Codon table of probabilities.
        """
        if nmm_frame_state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state` is NULL.")
        self._nmm_frame_state = nmm_frame_state
        self._baset = baset
        self._codont = codont
        alphabet = baset.alphabet
        super().__init__(lib.nmm_frame_state_super(self._nmm_frame_state), alphabet)

    @classmethod
    def create(
        cls: Type[FrameState],
        name: bytes,
        baset: BaseTable,
        codont: CodonTable,
        epsilon: float,
    ) -> FrameState:
        """
        Create frame state.

        Parameters
        ----------
        name
            State name.
        baset
            Base table of probabilities.
        codont
            Codon table of probabilities.
        epsilon
            Epsilon.
        """
        ptr = lib.nmm_frame_state_create(
            name, baset.nmm_base_table, codont.nmm_codon_table, epsilon
        )
        return FrameState(ptr, baset, codont)

    @property
    def base_table(self) -> BaseTable:
        return self._baset

    @property
    def codon_table(self) -> CodonTable:
        return self._codont

    def decode(self, seq: Sequence) -> Tuple[float, Codon]:
        state = self._nmm_frame_state
        any_symbol = self.alphabet.any_symbol
        codon = Codon.create(any_symbol * 3, self.alphabet)
        lprob = lib.nmm_frame_state_decode(state, seq.imm_seq, codon.nmm_codon)
        return lprob, codon

    @property
    def epsilon(self) -> float:
        return lib.nmm_frame_state_epsilon(self._nmm_frame_state)

    def __del__(self):
        if self._nmm_frame_state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._nmm_frame_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class CodonState(State[BaseAlphabet]):
    def __init__(self, nmm_codon_state: CData, codonp: CodonProb):
        """
        Codon state.

        Parameters
        ----------
        nmm_codon_state
            State pointer.
        codonp
            Codon probabilities.
        """
        if nmm_codon_state == ffi.NULL:
            raise RuntimeError("`nmm_codon_state` is NULL.")
        self._nmm_codon_state = nmm_codon_state
        self._codonp = codonp
        alphabet = codonp.alphabet
        super().__init__(lib.nmm_codon_state_super(nmm_codon_state), alphabet)

    @classmethod
    def create(cls: Type[CodonState], name: bytes, codonp: CodonProb) -> CodonState:
        """
        Create codon state.

        Parameters
        ----------
        name
            State name.
        codonp
            Codon probabilities.
        """
        ptr = lib.nmm_codon_state_create(name, codonp.nmm_codon_lprob)
        return CodonState(ptr, codonp)

    def __del__(self):
        if self._nmm_codon_state != ffi.NULL:
            lib.nmm_codon_state_destroy(self._nmm_codon_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
