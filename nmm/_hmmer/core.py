from typing import Sequence, Tuple

from .._hmm import HMM
from .._log import LOG1
from .._path import Path
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

            # hmm.del_state(core_nodes_trans[0][0].D)
            # hmm.del_state(core_nodes_trans[-1][0].I)

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
