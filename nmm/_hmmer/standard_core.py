from typing import Callable, List, Sequence


from .._hmm import HMM
from .._state import MuteState, NormalState
from .core import CoreModel, Node, NullModel, SpecialNode
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


class StandardCoreModel(CoreModel):
    def __init__(
        self, hmm: HMM, core_nodes: List[StandardNode], finalize: Callable[[], None]
    ):
        super().__init__(hmm, finalize)
        self._core_nodes = core_nodes

    def add_node(self, node: StandardNode, trans: Transitions):

        self._core_nodes.append(node)
        self._add_node(node, trans)

        if len(self._core_nodes) == 1:
            return

    def core_nodes(self) -> Sequence[StandardNode]:
        return self._core_nodes


class StandardNullModel(NullModel):
    def __init__(self, state: NormalState):
        super().__init__(state)
        self._normal_state = state

    @property
    def state(self) -> NormalState:
        return self._normal_state
