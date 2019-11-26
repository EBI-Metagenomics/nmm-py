from math import isnan
from ._path import Path, CPath
from ._log import LOG0
from ._state import State
from ._alphabet import Alphabet
from typing import Dict, NamedTuple


from ._ffi import ffi, lib

PathScore = NamedTuple("PathScore", [("score", float), ("path", CPath)])


class HMM:
    """
    Hidden Markov model.

    Parameters
    ----------
    alphabet : Alphabet
        Alphabet.
    """

    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._states: Dict[ffi.CData, State] = {}
        self._hmm = lib.imm_hmm_create(self._alphabet.imm_abc)
        if self._hmm == ffi.NULL:
            raise RuntimeError("`imm_hmm_create` failed.")

    def states(self) -> Dict[ffi.CData, State]:
        return self._states

    def set_start_lprob(self, state: State, lprob: float):
        err: int = lib.imm_hmm_set_start(self._hmm, state.imm_state, lprob)
        if err != 0:
            raise RuntimeError("Could not set start probability.")

    def transition(self, a: State, b: State):
        """
        Parameters
        ----------
        a : State
            Source state.
        b : State
            Destination state.
        """
        lprob: float = lib.imm_hmm_get_trans(self._hmm, a.imm_state, b.imm_state)
        if isnan(lprob):
            raise RuntimeError("Could not retrieve transition probability.")
        return lprob

    def set_transition(self, a: State, b: State, lprob: float):
        """
        Parameters
        ----------
        a : State
            Source state.
        b : State
            Destination state.
        lprob : float
            Transition probability in log-space.
        """
        if a.imm_state not in self._states:
            raise ValueError(f"State {a} not found.")

        if b.imm_state not in self._states:
            raise ValueError(f"State {b} not found.")

        err: int = lib.imm_hmm_set_trans(self._hmm, a.imm_state, b.imm_state, lprob)
        if err != 0:
            raise RuntimeError("Could not set transition probability.")

    @property
    def alphabet(self):
        return self._alphabet

    def add_state(self, state: State, start_lprob: float = LOG0):
        """
        Parameters
        ----------
        state
            Add state.
        start_lprob : bool
            Log-space probability of being the initial state.
        """
        err: int = lib.imm_hmm_add_state(self._hmm, state.imm_state, start_lprob)
        if err != 0:
            raise ValueError("Could not add state %s.", state)
        self._states[state.imm_state] = state

    def del_state(self, state: State):
        if state.imm_state not in self._states:
            raise ValueError(f"State {state} not found.")

        err: int = lib.imm_hmm_del_state(self._hmm, state.imm_state)
        if err != 0:
            raise RuntimeError(f"Could not delete state {state}.")

        del self._states[state.imm_state]

    def normalize(self):
        err: int = lib.imm_hmm_normalize(self._hmm)
        if err != 0:
            raise ValueError("Normalization error.")

    def normalize_transitions(self, state: State):
        err: int = lib.imm_hmm_normalize_trans(self._hmm, state.imm_state)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: bytes, path: Path):
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq, path.imm_path)
        if isnan(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(self, seq: bytes, end_state: State) -> PathScore:
        cpath = lib.imm_path_create()
        if cpath == ffi.NULL:
            raise RuntimeError("Could not create `cpath`.")
        try:
            lprob: float = lib.imm_hmm_viterbi(
                self._hmm, seq, end_state.imm_state, cpath
            )
        except Exception as e:
            lib.imm_path_destroy(cpath)
            raise e

        path = CPath(cpath)
        return PathScore(score=lprob, path=path)

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)
