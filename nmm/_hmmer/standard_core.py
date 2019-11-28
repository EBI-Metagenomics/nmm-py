from typing import Sequence, Tuple, Union, TypeVar, Type

from .._ffi import ffi, lib
from .._path import CPath, Path
from .._state import MuteState, NormalState
from .._step import Step
from .core import AltModel, Node, NullModel, SpecialNode
from .transition import Transitions


class StandardNode(Node):
    def __init__(self, M: NormalState, I: NormalState, D: MuteState):
        self._M = M
        self._I = I
        self._D = D

    @property
    def M(self) -> NormalState:
        return self._M

    @property
    def I(self) -> NormalState:
        return self._I

    @property
    def D(self) -> MuteState:
        return self._D


class StandardSpecialNode(SpecialNode):
    def __init__(
        self,
        S: MuteState,
        N: NormalState,
        B: MuteState,
        E: MuteState,
        J: NormalState,
        C: NormalState,
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
    def N(self) -> NormalState:
        return self._N

    @property
    def B(self) -> MuteState:
        return self._B

    @property
    def E(self) -> MuteState:
        return self._E

    @property
    def J(self) -> NormalState:
        return self._J

    @property
    def C(self) -> NormalState:
        return self._C

    @property
    def T(self) -> MuteState:
        return self._T


class StandardNullModel(NullModel):
    def __init__(self, state: NormalState):
        super().__init__(state)
        self._normal_state = state

    @property
    def state(self) -> NormalState:
        return self._normal_state


class StandardAltModel(AltModel):
    def __init__(
        self,
        special_node: StandardSpecialNode,
        nodes_trans: Sequence[Tuple[StandardNode, Transitions]],
    ):
        self._special_node = special_node
        self._core_nodes = [nt[0] for nt in nodes_trans]
        super().__init__(special_node, nodes_trans)

    @property
    def length(self):
        return len(self._core_nodes)

    def core_nodes(self) -> Sequence[StandardNode]:
        return self._core_nodes

    @property
    def special_node(self) -> StandardSpecialNode:
        return self._special_node

    def viterbi(self, seq: bytes) -> Tuple[float, CPath]:
        return self._hmm.viterbi(seq, self.special_node.T)


class StandardStep(Step):
    def __init__(
        self, imm_step: ffi.CData, state: Union[MuteState, NormalState], seq_len: int
    ):
        super().__init__(imm_step, state, seq_len)
        self._state = state

    @property
    def state(self) -> Union[MuteState, NormalState]:
        return self._state


# T = TypeVar("T", bound="StandardPath")


# class StandardPath(Path):
#     def __init__(self):
#         super().__init__()

#     @classmethod
#     def create(
#         cls: Type[T], steps: Sequence[Tuple[Union[MuteState, NormalState], int]]
#     ) -> T:
#         path = cls()
#         for state, seq_len in steps:
#             path.append(state, seq_len)
#         return path

# def append(self, state: Union[MuteState, NormalState], seq_len: int) -> ffi.CData:
#     imm_step = super().append(state.imm_state, seq_len)
#     step = Step(imm_step, state, seq_len)
#     self.__steps.append(step)
#     return step
