from math import isnan
from ._string import make_sure_bytes
from ._path import Path
from ._log import LOG
from ._state import State
from ._alphabet import Alphabet
from typing import Dict

from ._ffi import ffi, lib


class HMM:
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._hmm = lib.imm_hmm_create(self._alphabet.cdata)
        self._states: Dict[ffi.CData, State] = {}

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)

    @property
    def states(self):
        return self._states

    def set_start_lprob(self, state: State, lprob: float):
        err: int = lib.imm_hmm_set_start_lprob(self._hmm, state.cdata, lprob)
        if err != 0:
            raise ValueError("Could not set start probability.")

    def trans(self, a: State, b: State):
        """
        Parameters
        ----------
        a : State
            Source state.
        b : State
            Destination state.
        """
        lprob: float = lib.imm_hmm_get_trans(self._hmm, a.cdata, b.cdata)
        if isnan(lprob):
            raise ValueError("Could not retrieve transition probability.")
        return lprob

    def set_trans(self, a: State, b: State, lprob: float):
        """
        Parameters
        ----------
        a : State
            Source state name.
        b : State
            Destination state name.
        lprob : float
            Transition probability in log-space.
        """
        err: int = lib.imm_hmm_set_trans(self._hmm, a.cdata, b.cdata, lprob)
        if err != 0:
            raise ValueError("Could not set transition probability.")

    @property
    def alphabet(self):
        return self._alphabet

    def add_state(self, state: State, start_lprob: float = LOG(0.0)):
        """
        Parameters
        ----------
        state
            Add state.
        start_lprob : bool
            Log-space probability of being the initial state.
        """
        err: int = lib.imm_hmm_add_state(self._hmm, state.cdata, start_lprob)
        if err != 0:
            raise ValueError("Could not add state %s.", state)
        self._states[state.cdata] = state

    def del_state(self, state: State):
        if state.cdata not in self._states:
            raise ValueError(f"State {state} not found.")

        err: int = lib.imm_hmm_del_state(self._hmm, state.cdata)
        if err != 0:
            raise ValueError(f"Could not delete state {state}.")

        del self._states[state.cdata]

    def normalize(self):
        err: int = lib.imm_hmm_normalize(self._hmm)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: str, path: Path):
        seq = make_sure_bytes(seq)
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq, path.cdata)
        if isnan(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(self, seq: str, end_state: State):
        seq = make_sure_bytes(seq)
        lprob: float = lib.imm_hmm_viterbi(self._hmm, seq, end_state.cdata)
        if isnan(lprob):
            raise ValueError("Could not calculate the viterbi score.")
        return lprob
