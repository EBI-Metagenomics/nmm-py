from ._sequence import CSequence
from ._path import CPath
from ._state import CState
from ._alphabet import CAlphabet
from ._lprob import LPROB_ZERO, lprob_is_valid
from ._results import CResults
from typing import Dict, Tuple


from ._ffi import ffi, lib


class HMM:
    """
    Hidden Markov model.

    Parameters
    ----------
    alphabet : `CAlphabet`
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
        if lib.imm_hmm_set_start(self._hmm, state.imm_state, lprob) != 0:
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
        if not lprob_is_valid(lprob):
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

    def add_state(self, state: CState, start_lprob: float = LPROB_ZERO):
        """
        Parameters
        ----------
        state :  `CState`
            Add state.
        start_lprob : float
            Log-space probability of being the initial state.
        """
        if lib.imm_hmm_add_state(self._hmm, state.imm_state, start_lprob) != 0:
            raise ValueError(f"Could not add state {str(state.name)}.")
        self._states[state.imm_state] = state

    def del_state(self, state: CState):
        if state.imm_state not in self._states:
            raise ValueError(f"State {state} not found.")

        err: int = lib.imm_hmm_del_state(self._hmm, state.imm_state)
        if err != 0:
            raise RuntimeError(f"Could not delete state {state}.")

        del self._states[state.imm_state]

    def normalize(self):
        if lib.imm_hmm_normalize(self._hmm) != 0:
            raise ValueError("Normalization error.")

    def normalize_transitions(self, state: CState):
        err: int = lib.imm_hmm_normalize_trans(self._hmm, state.imm_state)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: CSequence, path: CPath):
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq.imm_seq, path.imm_path)
        if not lprob_is_valid(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(self, seq: CSequence, end_state: CState, window_length: int = 0):
        state = end_state.imm_state
        r = lib.imm_hmm_viterbi(self._hmm, seq.imm_seq, state, window_length)
        return CResults(r, seq)

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)
