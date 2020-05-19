from __future__ import annotations
from typing import Tuple, Type, TypeVar, Dict

from enum import Enum
from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import (
    Sequence,
    State,
    wrap_imm_state as imm_wrap_imm_state,
)
from ._codon_prob import CodonProb
from ._base_alphabet import BaseAlphabet
from ._base_table import BaseTable
from ._codon import Codon
from ._codon_table import CodonTable
from .._imm import Alphabet

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


def wrap_imm_state(
    imm_state: CData,
    alphabet: T,
    base_tables: Dict[CData, BaseTable],
    codon_tables: Dict[CData, CodonTable],
    codon_probs: Dict[CData, CodonProb],
) -> State:
    try:
        state_type = StateType(lib.imm_state_type_id(imm_state))
    except ValueError:
        return imm_wrap_imm_state(imm_state, alphabet)
    if state_type == StateType.CODON:
        pass
    elif state_type == StateType.FRAME:
        nmm_frame_state = lib.nmm_frame_state_derived(imm_state)
        if nmm_frame_state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state` is NULL.")

        nmm_base_table = lib.nmm_frame_state_baset(nmm_frame_state)
        nmm_codon_table = lib.nmm_frame_state_codont(nmm_frame_state)

        baset = base_tables[nmm_base_table]
        codont = codon_tables[nmm_codon_table]

        return FrameState(nmm_frame_state, baset, codont)
    raise ValueError(f"Unknown state type: {lib.imm_state_type_id(imm_state)}.")
