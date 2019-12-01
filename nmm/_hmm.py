from math import isnan
from ._path import CPath
from ._log import LOG0
from ._state import CState
from ._alphabet import CAlphabet
from typing import Dict, Tuple


from ._ffi import ffi, lib


class HMM:
    """
    Hidden Markov model.

    Parameters
    ----------
    alphabet : CAlphabet
        Alphabet.
    """

    def __init__(self, alphabet: CAlphabet):
        self._alphabet = alphabet
        self._states: Dict[ffi.CData, CState] = {}
        self._hmm = lib.imm_hmm_create(self._alphabet.imm_abc)
        if self._hmm == ffi.NULL:
            raise RuntimeError("`imm_hmm_create` failed.")

    def states(self) -> Dict[ffi.CData, CState]:
        return self._states

    def set_start_lprob(self, state: CState, lprob: float):
        err: int = lib.imm_hmm_set_start(self._hmm, state.imm_state, lprob)
        if err != 0:
            raise RuntimeError("Could not set start probability.")

    def transition(self, a: CState, b: CState):
        """
        Parameters
        ----------
        a : CState
            Source state.
        b : CState
            Destination state.
        """
        lprob: float = lib.imm_hmm_get_trans(self._hmm, a.imm_state, b.imm_state)
        if isnan(lprob):
            raise RuntimeError("Could not retrieve transition probability.")
        return lprob

    def set_transition(self, a: CState, b: CState, lprob: float):
        """
        Parameters
        ----------
        a : CState
            Source state.
        b : CState
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
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    def add_state(self, state: CState, start_lprob: float = LOG0):
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

    def del_state(self, state: CState):
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

    def normalize_transitions(self, state: CState):
        err: int = lib.imm_hmm_normalize_trans(self._hmm, state.imm_state)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: bytes, path: CPath):
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq, path.imm_path)
        if isnan(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(self, seq: bytes, end_state: CState) -> Tuple[float, CPath]:
        imm_path = lib.imm_path_create()
        if imm_path == ffi.NULL:
            raise RuntimeError("Could not create `imm_path`.")
        try:
            imm_state = end_state.imm_state
            lprob: float = lib.imm_hmm_viterbi(self._hmm, seq, imm_state, imm_path)
        except Exception as e:
            lib.imm_path_destroy(imm_path)
            raise e

        path = CPath.create_cpath()
        for step in CPath(imm_path).steps():
            path.append_cstep(self._states[step.state.imm_state], step.seq_len)

        return (lprob, path)

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)
