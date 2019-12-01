from typing import Dict, Iterator, List, Sequence, Tuple, TypeVar, Union

from .._ffi import ffi
from .._path import CPath
from .._state import MuteState, FrameState
from .._step import CStep
from .core import AltModel, Node, NullModel, SpecialNode
from .transition import Transitions


class FrameStep(CStep):
    def __init__(self, imm_step: ffi.CData, state: Union[MuteState, FrameState]):
        super().__init__(imm_step)
        self._state = state

    @property
    def state(self) -> Union[MuteState, FrameState]:
        return self._state


T = TypeVar("T", bound="FramePath")


class FramePath(CPath):
    def __init__(self):
        super().__init__()
        self._steps: List[FrameStep] = []

    def append_frame_step(
        self, state: Union[MuteState, FrameState], seq_len: int
    ) -> FrameStep:
        cstep = self.append_cstep(state, seq_len)
        self._steps.append(FrameStep(cstep.imm_step, state))
        return self._steps[-1]

    def steps(self) -> Iterator[FrameStep]:
        return iter(self._steps)


class FrameNode(Node):
    def __init__(self, M: FrameState, I: FrameState, D: MuteState):
        self._M = M
        self._I = I
        self._D = D

    @property
    def M(self) -> FrameState:
        return self._M

    @property
    def I(self) -> FrameState:
        return self._I

    @property
    def D(self) -> MuteState:
        return self._D

    def states(self) -> List[Union[MuteState, FrameState]]:
        return [self._M, self._I, self._D]


class FrameSpecialNode(SpecialNode):
    def __init__(
        self,
        S: MuteState,
        N: FrameState,
        B: MuteState,
        E: MuteState,
        J: FrameState,
        C: FrameState,
        T: MuteState,
    ):
        self._S = S
        self._N = N
        self._B = B
        self._E = E
        self._J = J
        self._C = C
        self._T = T

    @property
    def S(self) -> MuteState:
        return self._S

    @property
    def N(self) -> FrameState:
        return self._N

    @property
    def B(self) -> MuteState:
        return self._B

    @property
    def E(self) -> MuteState:
        return self._E

    @property
    def J(self) -> FrameState:
        return self._J

    @property
    def C(self) -> FrameState:
        return self._C

    @property
    def T(self) -> MuteState:
        return self._T

    def states(self) -> List[Union[MuteState, FrameState]]:
        return [self._S, self._N, self._B, self._E, self._J, self._C, self._T]


class FrameNullModel(NullModel):
    def __init__(self, state: FrameState):
        super().__init__(state)
        self._normal_state = state

    @property
    def state(self) -> FrameState:
        return self._normal_state


class FrameAltModel(AltModel):
    def __init__(
        self,
        special_node: FrameSpecialNode,
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        self._special_node = special_node
        self._core_nodes = [nt[0] for nt in nodes_trans]
        self._states: Dict[ffi.CData, Union[MuteState, FrameState]] = {}

        for node in self._core_nodes:
            for state in node.states():
                self._states[state.imm_state] = state

        for state in special_node.states():
            self._states[state.imm_state] = state

        super().__init__(special_node, nodes_trans)

    @property
    def length(self):
        return len(self._core_nodes)

    def core_nodes(self) -> Sequence[FrameNode]:
        return self._core_nodes

    @property
    def special_node(self) -> FrameSpecialNode:
        return self._special_node

    def viterbi(self, seq: bytes) -> Tuple[float, FramePath]:
        score, path = super().viterbi(seq)

        spath = FramePath()
        for step in path.steps():
            imm_state = step.state.imm_state
            spath.append_frame_step(self._states[imm_state], step.seq_len)

        return (score, spath)
