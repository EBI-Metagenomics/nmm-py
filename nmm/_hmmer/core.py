from typing import Callable, List, NamedTuple

from .._hmm import HMM
from .._log import LOG1
from .._path import Path
from .._state import State
from .transition import Transitions

Node = NamedTuple("Node", [("M", State), ("I", State), ("D", State)])


class CoreModel:
    def __init__(self, hmm: HMM, core_nodes: List[Node], finalize: Callable[[], None]):
        self._hmm = hmm
        self._core_nodes = core_nodes
        self._finalize = finalize

    def add_node(self, node: Node, trans: Transitions):
        self._hmm.add_state(node.M)
        self._hmm.add_state(node.I)
        self._hmm.add_state(node.D)

        self._core_nodes.append(node)

        if len(self._core_nodes) == 1:
            return

        prev = self._core_nodes[-2]
        self._hmm.set_transition(prev.M, node.M, trans.MM)
        self._hmm.set_transition(prev.M, prev.I, trans.MI)
        self._hmm.set_transition(prev.M, node.D, trans.MD)
        self._hmm.set_transition(prev.I, node.M, trans.IM)
        self._hmm.set_transition(prev.I, prev.I, trans.II)
        self._hmm.set_transition(prev.D, node.M, trans.DM)
        self._hmm.set_transition(prev.D, node.D, trans.DD)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if len(self._core_nodes) > 0:
            self._hmm.del_state(self._core_nodes[0].D)
            self._hmm.del_state(self._core_nodes[-1].I)
        self._finalize()
        del type
        del value
        del traceback


class NullModel:
    def __init__(self, state: State):
        self._state = state
        self._hmm = HMM(state.alphabet)
        self._hmm.add_state(state, LOG1)

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def set_transition(self, lprob: float):
        self._hmm.set_transition(self._state, self._state, lprob)

    def likelihood(self, seq: bytes):
        path = Path.create([(self._state, 1) for i in range(len(seq))])
        return self._hmm.likelihood(seq, path)
