from __future__ import annotations

from typing import Dict, Generic, Optional, Type, TypeVar

from .._cdata import CData
from .._ffi import ffi, lib
from ._alphabet import Alphabet
from ._lprob import lprob_is_valid, lprob_is_zero, lprob_zero
from ._path import Path
from ._sequence import Sequence
from ._state import State

TState = TypeVar("TState", bound=State)


class HMM(Generic[TState]):
    """
    Hidden Markov model.

    Parameters
    ----------
    alphabet : `Alphabet`
        Alphabet.
    """

    def __init__(
        self,
        imm_hmm: CData,
        alphabet: Alphabet,
        states: Optional[Dict[CData, TState]] = None,
    ):
        if imm_hmm == ffi.NULL:
            raise RuntimeError("`imm_hmm` is NULL.")
        self._alphabet = alphabet
        if states is None:
            self._states: Dict[CData, TState] = {}
        else:
            self._states = states
        self._imm_hmm = imm_hmm

    @classmethod
    def create(cls: Type[HMM], alphabet: Alphabet) -> HMM:
        """
        Create HMM.

        Parameters
        ----------
        alphabet
            Alphabet.
        """
        imm_hmm = lib.imm_hmm_create(alphabet.imm_abc)
        return cls(imm_hmm, alphabet)

    @property
    def imm_hmm(self) -> CData:
        return self._imm_hmm

    def states(self) -> Dict[CData, TState]:
        return self._states

    def set_start_lprob(self, state: TState, lprob: float):
        if lib.imm_hmm_set_start(self._imm_hmm, state.imm_state, lprob) != 0:
            raise RuntimeError("Could not set start probability.")

    def transition(self, a: TState, b: TState):
        """
        Parameters
        ----------
        a
            Source state.
        b
            Destination state.
        """
        lprob: float = lib.imm_hmm_get_trans(self._imm_hmm, a.imm_state, b.imm_state)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not retrieve transition probability.")
        return lprob

    def set_transition(self, a: TState, b: TState, lprob: float):
        """
        Parameters
        ----------
        a
            Source state.
        b
            Destination state.
        lprob : float
            Transition probability in log-space.
        """
        if a.imm_state not in self._states:
            raise ValueError(f"State {a} not found.")

        if b.imm_state not in self._states:
            raise ValueError(f"State {b} not found.")

        err: int = lib.imm_hmm_set_trans(self._imm_hmm, a.imm_state, b.imm_state, lprob)
        if err != 0:
            raise RuntimeError("Could not set transition probability.")

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    def add_state(self, state: TState, start_lprob: float = lprob_zero()):
        """
        Parameters
        ----------
        state
            Add state.
        start_lprob
            Log-space probability of being the initial state.
        """
        if lib.imm_hmm_add_state(self._imm_hmm, state.imm_state, start_lprob) != 0:
            raise ValueError(f"Could not add state {str(state.name)}.")
        self._states[state.imm_state] = state

    def del_state(self, state: TState):
        if state.imm_state not in self._states:
            raise ValueError(f"State {state} not found.")

        err: int = lib.imm_hmm_del_state(self._imm_hmm, state.imm_state)
        if err != 0:
            raise RuntimeError(f"Could not delete state {state}.")

        del self._states[state.imm_state]

    def normalize(self):
        if lib.imm_hmm_normalize(self._imm_hmm) != 0:
            raise ValueError("Normalization error.")

    def normalize_transitions(self, state: TState):
        err: int = lib.imm_hmm_normalize_trans(self._imm_hmm, state.imm_state)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: Sequence, path: Path) -> float:
        lprob: float = lib.imm_hmm_likelihood(self._imm_hmm, seq.imm_seq, path.imm_path)
        if not lprob_is_valid(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def create_dp(self, end_state: TState):
        from ._dp import DP

        imm_state = end_state.imm_state
        imm_dp = lib.imm_hmm_create_dp(self._imm_hmm, imm_state)
        if imm_dp == ffi.NULL:
            raise RuntimeError("Could not create dp.")

        return DP(imm_dp, self)

    def view(self):
        from graphviz import Digraph

        dot = Digraph(comment="HMM")

        for state in self._states.values():
            dot.node(state.name.decode(), state.name.decode())

        for state0 in self._states.values():
            for state1 in self._states.values():
                t = self.transition(state0, state1)
                if lprob_is_zero(t):
                    continue
                label = f"{t:.6f}"
                dot.edge(state0.name.decode(), state1.name.decode(), label=label)

        dot.view()

    def __del__(self):
        if self._imm_hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._imm_hmm)
