from ._sequence import Sequence
from ._path import Path
from ._state import CState
from ._alphabet import Alphabet
from ._lprob import LPROB_ZERO, lprob_is_valid
from ._results import CResults
from typing import Dict


from .._ffi import ffi, lib


class HMM:
    """
    Hidden Markov model.

    Parameters
    ----------
    alphabet : `Alphabet`
        Alphabet.
    """

    def __init__(self, alphabet: Alphabet):
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
    def alphabet(self) -> Alphabet:
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

    def likelihood(self, seq: Sequence, path: Path):
        lprob: float = lib.imm_hmm_likelihood(self._hmm, seq.imm_seq, path.imm_path)
        if not lprob_is_valid(lprob):
            raise ValueError("Could not calculate the likelihood.")
        return lprob

    def viterbi(
        self, seq: Sequence, end_state: CState, window_length: int = 0
    ) -> CResults:
        from ._results import wrap_imm_results

        imm_seq = seq.imm_seq
        imm_state = end_state.imm_state

        imm_results = lib.imm_hmm_viterbi(self._hmm, imm_seq, imm_state, window_length)
        if imm_results == ffi.NULL:
            raise RuntimeError("Could not run viterbi.")

        return wrap_imm_results(imm_results, seq, self._states)

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)
