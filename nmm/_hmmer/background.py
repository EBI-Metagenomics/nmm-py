from .._path import Path
from .._hmm import HMM
from .._state import NormalState


class BackgroundModel:
    def __init__(self, state: NormalState):
        self._state = state
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(self._state, 0.0)

    @property
    def state(self):
        return self._state

    @property
    def hmm(self):
        return self._hmm

    def set_trans(self, lprob: float):
        self._hmm.set_trans(self._state, self._state, lprob)

    def likelihood(self, seq: str):
        s = seq.encode()
        path = [(self._state, 1)] * len(s)
        return self._hmm.likelihood(seq, Path(path))
