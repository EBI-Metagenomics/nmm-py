from typing import Callable, Sequence, Tuple

from math import log
from .._hmm import HMM
from .._log import LOG1, LOG0
from .._path import Path, CPath
from .._state import State, MuteState
from .transition import Transitions, SpecialTransitions


class Node:
    @property
    def M(self) -> State:
        raise NotImplementedError()

    @property
    def I(self) -> State:
        raise NotImplementedError()

    @property
    def D(self) -> State:
        raise NotImplementedError()


class SpecialNode:
    @property
    def S(self) -> MuteState:
        raise NotImplementedError()

    @property
    def N(self) -> State:
        raise NotImplementedError()

    @property
    def B(self) -> MuteState:
        raise NotImplementedError()

    @property
    def E(self) -> MuteState:
        raise NotImplementedError()

    @property
    def J(self) -> State:
        raise NotImplementedError()

    @property
    def C(self) -> State:
        raise NotImplementedError()

    @property
    def T(self) -> MuteState:
        raise NotImplementedError()


class CoreModel:
    def __init__(self, hmm: HMM, finalize: Callable[[], None]):
        self._hmm = hmm
        self._finalize = finalize

    def _add_node(self, node: Node, trans: Transitions):
        self._hmm.add_state(node.M)
        self._hmm.add_state(node.I)
        self._hmm.add_state(node.D)

        if self.model_length == 1:
            return

        prev = self.core_nodes()[-2]
        self._hmm.set_transition(prev.M, node.M, trans.MM)
        self._hmm.set_transition(prev.M, prev.I, trans.MI)
        self._hmm.set_transition(prev.M, node.D, trans.MD)
        self._hmm.set_transition(prev.I, node.M, trans.IM)
        self._hmm.set_transition(prev.I, prev.I, trans.II)
        self._hmm.set_transition(prev.D, node.M, trans.DM)
        self._hmm.set_transition(prev.D, node.D, trans.DD)

    def core_nodes(self) -> Sequence[Node]:
        raise NotImplementedError()

    @property
    def model_length(self) -> int:
        return len(self.core_nodes())

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.model_length > 0:
            self._hmm.del_state(self.core_nodes()[0].D)
            self._hmm.del_state(self.core_nodes()[-1].I)
        self._finalize()
        del type
        del value
        del traceback


class NullModel:
    def __init__(self, state: State):
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(state, LOG1)

    @property
    def state(self) -> State:
        raise NotImplementedError()

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def set_transition(self, lprob: float):
        self._hmm.set_transition(self.state, self.state, lprob)

    def likelihood(self, seq: bytes):
        path = Path.create([(self.state, 1) for i in range(len(seq))])
        return self._hmm.likelihood(seq, path)


class Profile:
    def __init__(self, special_node: SpecialNode):

        hmm = HMM(special_node.S.alphabet)
        hmm.add_state(special_node.S, LOG1)
        hmm.add_state(special_node.N)
        hmm.add_state(special_node.B)
        hmm.add_state(special_node.E)
        hmm.add_state(special_node.J)
        hmm.add_state(special_node.C)
        hmm.add_state(special_node.T)

        self._hmm = hmm
        self._multiple_hits: bool = True
        self._special_transitions = SpecialTransitions()

    @property
    def null_model(self) -> NullModel:
        raise NotImplementedError()

    @property
    def special_node(self) -> SpecialNode:
        raise NotImplementedError()

    @property
    def length(self) -> int:
        raise NotImplementedError()

    def core_nodes(self) -> Sequence[Node]:
        raise NotImplementedError()

    @property
    def core_model(self):
        raise NotImplementedError()

    def set_multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def _finalize(self):
        self._set_fragment_length()

    def _set_fragment_length(self):
        if self.length == 0:
            return

        B = self._special_node.B
        E = self._special_node.E

        # Uniform local alignment fragment length distribution
        t = self._special_transitions
        t.BM = log(2) - log(self.length) - log(self.length + 1)
        t.ME = 0.0
        for node in self.core_nodes():
            self._hmm.set_transition(B, node.M, t.BM)
            self._hmm.set_transition(node.M, E, t.ME)

        for node in self.core_nodes()[1:]:
            self._hmm.set_transition(node.D, E, 0.0)

    def _set_target_length(self, seq: bytes):
        from math import exp

        L = len(seq)
        if L == 0:
            return

        if self._multiple_hits:
            lq = -log(2)
        else:
            lq = LOG0

        q = exp(lq)
        lp = log(L) - log(L + 2 + q / (1 - q))
        l1p = log(2 + q / (1 - q)) - log(L + 2 + q / (1 - q))
        lr = log(L) - log(L + 1)

        t = self._special_transitions

        t.NN = t.CC = t.JJ = lp
        t.NB = t.CT = t.JB = l1p
        t.RR = lr
        t.EC = t.EJ = lq

        node = self.special_node

        self._hmm.set_transition(node.S, node.B, t.NB)
        self._hmm.set_transition(node.S, node.N, t.NN)
        self._hmm.set_transition(node.N, node.N, t.NN)
        self._hmm.set_transition(node.N, node.B, t.NB)

        self._hmm.set_transition(node.E, node.T, t.EC + t.CT)
        self._hmm.set_transition(node.E, node.C, t.EC + t.CC)
        self._hmm.set_transition(node.C, node.C, t.CC)
        self._hmm.set_transition(node.C, node.T, t.CT)

        self._hmm.set_transition(node.E, node.B, t.EJ + t.JB)
        self._hmm.set_transition(node.E, node.J, t.EJ + t.JJ)
        self._hmm.set_transition(node.J, node.J, t.JJ)
        self._hmm.set_transition(node.J, node.B, t.JB)

        self.null_model.set_transition(t.RR)

    def _viterbi(self, seq: bytes) -> Tuple[float, CPath]:
        self._set_target_length(seq)
        return self._hmm.viterbi(seq, self.special_node.T)
