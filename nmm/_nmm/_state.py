from __future__ import annotations
from typing import Dict, Tuple, Type

from enum import Enum
from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import (
    Alphabet,
    Sequence,
    SequenceTable,
    State,
    TableState,
    wrap_imm_state as imm_wrap_imm_state,
)
from ._codon_prob import CodonProb
from ._base_alphabet import BaseAlphabet
from ._base_table import BaseTable
from ._codon import Codon
from ._codon_table import CodonTable


class StateType(Enum):
    CODON = 0x10
    FRAME = 0x11


class FrameState(State[BaseAlphabet]):
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
            name, baset.nmm_base_table, codont.nmm_codon_table, epsilon
        )
        if state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state_create` failed.")

        self._baset = baset
        self._codont = codont
        self._epsilon = epsilon
        self._nmm_frame_state = state
        alphabet = baset.alphabet
        super().__init__(lib.nmm_frame_state_super(self._nmm_frame_state), alphabet)

    def decode(self, seq: Sequence) -> Tuple[float, Codon]:
        state = self._nmm_frame_state
        any_symbol = self.alphabet.any_symbol
        codon = Codon.create(any_symbol * 3, self.alphabet)
        lprob = lib.nmm_frame_state_decode(state, seq.imm_seq, codon.nmm_codon)
        return lprob, codon

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


# states: Mapping[CData, TState]
def wrap_imm_state(imm_state: CData) -> State:
    state_type: int = lib.imm_state_type_id(imm_state)
    if state_type == StateType.CODON:
        pass
    elif state_type == StateType.FRAME:
        pass
    return imm_wrap_imm_state(imm_state)
