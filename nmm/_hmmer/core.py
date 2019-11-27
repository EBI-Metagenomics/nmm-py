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


class NullModel:
    def __init__(self, state: State):
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(state, LOG1)

    @property
    def state(self) -> State:
        raise NotImplementedError()

    # @property
    # def hmm(self) -> HMM:
    #     return self._hmm

    def set_transition(self, lprob: float):
        self._hmm.set_transition(self.state, self.state, lprob)

    def likelihood(self, seq: bytes):
        path = Path.create([(self.state, 1) for i in range(len(seq))])
        return self._hmm.likelihood(seq, path)


class AltModel:
    def __init__(
        self,
        special_node: SpecialNode,
        core_nodes_trans: Sequence[Tuple[Node, Transitions]],
    ):
        hmm = HMM(special_node.S.alphabet)
        hmm.add_state(special_node.S, LOG1)
        hmm.add_state(special_node.N)
        hmm.add_state(special_node.B)
        hmm.add_state(special_node.E)
        hmm.add_state(special_node.J)
        hmm.add_state(special_node.C)
        hmm.add_state(special_node.T)

        self._special_transitions = SpecialTransitions()

        if len(core_nodes_trans) > 0:
            node, trans = core_nodes_trans[0]
            hmm.add_state(node.M)
            hmm.add_state(node.I)
            hmm.add_state(node.D)
            prev = node

            for node, trans in core_nodes_trans[1:]:
                hmm.add_state(node.M)
                hmm.add_state(node.I)
                hmm.add_state(node.D)

                hmm.set_transition(prev.M, node.M, trans.MM)
                hmm.set_transition(prev.M, prev.I, trans.MI)
                hmm.set_transition(prev.M, node.D, trans.MD)
                hmm.set_transition(prev.I, node.M, trans.IM)
                hmm.set_transition(prev.I, prev.I, trans.II)
                hmm.set_transition(prev.D, node.M, trans.DM)
                hmm.set_transition(prev.D, node.D, trans.DD)
                prev = node

            hmm.del_state(core_nodes_trans[0][0].D)
            hmm.del_state(core_nodes_trans[-1][0].I)

        self._hmm = hmm

    def set_transition(self, a: State, b: State, lprob: float):
        self._hmm.set_transition(a, b, lprob)

    def core_nodes(self) -> Sequence[Node]:
        raise NotImplementedError()

    @property
    def special_node(self) -> SpecialNode:
        raise NotImplementedError()

    @property
    def special_transitions(self) -> SpecialTransitions:
        return self._special_transitions

    @property
    def length(self) -> int:
        raise NotImplementedError()


class Profile:
    def __init__(self):
        self._multiple_hits: bool = True

    @property
    def null_model(self) -> NullModel:
        raise NotImplementedError()

    @property
    def alt_model(self) -> AltModel:
        raise NotImplementedError()

    @property
    def multiple_hits(self) -> bool:
        return self._multiple_hits

    @multiple_hits.setter
    def multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    # @property
    # def hmm(self) -> HMM:
    #     return self._hmm

    # def _finalize(self):
    #     self._set_fragment_length()

    def _set_fragment_length(self):
        if self.alt_model.length == 0:
            return

        B = self.alt_model.special_node.B
        E = self.alt_model.special_node.E

        # Uniform local alignment fragment length distribution
        t = self.alt_model.special_transitions
        t.BM = log(2) - log(self.alt_model.length) - log(self.alt_model.length + 1)
        t.ME = 0.0
        for node in self.alt_model.core_nodes():
            self.alt_model.set_transition(B, node.M, t.BM)
            self.alt_model.set_transition(node.M, E, t.ME)

        for node in self.alt_model.core_nodes()[1:]:
            self.alt_model.set_transition(node.D, E, 0.0)

    def _set_target_length(self, length: int):
        from math import exp

        L = length
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

        t = self.alt_model.special_transitions

        t.NN = t.CC = t.JJ = lp
        t.NB = t.CT = t.JB = l1p
        t.RR = lr
        t.EC = t.EJ = lq

        node = self.alt_model.special_node

        self.alt_model.set_transition(node.S, node.B, t.NB)
        self.alt_model.set_transition(node.S, node.N, t.NN)
        self.alt_model.set_transition(node.N, node.N, t.NN)
        self.alt_model.set_transition(node.N, node.B, t.NB)

        self.alt_model.set_transition(node.E, node.T, t.EC + t.CT)
        self.alt_model.set_transition(node.E, node.C, t.EC + t.CC)
        self.alt_model.set_transition(node.C, node.C, t.CC)
        self.alt_model.set_transition(node.C, node.T, t.CT)

        self.alt_model.set_transition(node.E, node.B, t.EJ + t.JB)
        self.alt_model.set_transition(node.E, node.J, t.EJ + t.JJ)
        self.alt_model.set_transition(node.J, node.J, t.JJ)
        self.alt_model.set_transition(node.J, node.B, t.JB)

        self.null_model.set_transition(t.RR)
