from math import isnan
from ._path import Path
from ._log import LOG0
from ._state import State
from ._alphabet import Alphabet
from bidict import bidict
from typing import Dict, Optional

from ._ffi import ffi, lib


class HMM:
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._states: Dict[ffi.CData, State] = {}
        self._state_names = bidict()
        self._hmm = lib.imm_hmm_create(self._alphabet.cdata)

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)

    @property
    def states(self):
        return self._states

    def set_start_lprob(self, state: State, lprob: float):
        err: int = lib.imm_hmm_set_start_lprob(self._hmm, state.cdata, lprob)
        if err != 0:
            raise RuntimeError("Could not set start probability.")

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
            raise RuntimeError("Could not retrieve transition probability.")
        return lprob

    def set_trans(self, a: State, b: State, lprob: float):
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
        if a.cdata not in self._states:
            raise ValueError(f"State {a} not found.")

        if b.cdata not in self._states:
            raise ValueError(f"State {b} not found.")

        err: int = lib.imm_hmm_set_trans(self._hmm, a.cdata, b.cdata, lprob)
        if err != 0:
            raise RuntimeError("Could not set transition probability.")

    def find_state(self, state_name: str):
        return self._states[self._state_names.inverse[state_name]]

    @property
    def alphabet(self):
        return self._alphabet

    def add_state(
        self, state: State, start_lprob: float = LOG0, name: Optional[str] = None
    ):
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
        if name is not None:
            self._state_names[state.cdata] = name

    def del_state(self, state: State):
        if state.cdata not in self._states:
            raise ValueError(f"State {state} not found.")

        err: int = lib.imm_hmm_del_state(self._hmm, state.cdata)
        if err != 0:
            raise RuntimeError(f"Could not delete state {state}.")

        del self._states[state.cdata]
        try:
            self._state_names.pop(state.cdata)
        except KeyError:
            pass

    def normalize(self):
        err: int = lib.imm_hmm_normalize(self._hmm)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: str, path: Path):
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq.encode(), path.cdata)
        if isnan(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(self, seq: str, end_state: State):
        lprob: float = lib.imm_hmm_viterbi(self._hmm, seq.encode(), end_state.cdata)
        if isnan(lprob):
            raise ValueError("Could not calculate the viterbi score.")
        return lprob
