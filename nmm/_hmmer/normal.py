from math import log
from typing import List, NamedTuple

from hmmer_reader import HMMEReader

from .._alphabet import Alphabet
from .._hmm import HMM, PathScore
from .._log import LOG0, LOG1
from .._state import MuteState, NormalState
from .core import CoreModel, NullModel, SpecialTrans, Trans
from .path import HMMERResult

Node = NamedTuple("Node", [("M", NormalState), ("I", NormalState), ("D", MuteState)])

SpecialNode = NamedTuple(
    "SpecialNode",
    [
        ("S", MuteState),
        ("N", NormalState),
        ("B", MuteState),
        ("E", MuteState),
        ("J", NormalState),
        ("C", NormalState),
        ("T", MuteState),
    ],
)


class NormalNullModel(NullModel):
    def __init__(self, state: NormalState):
        super().__init__(state)
        self._normal_state = state

    @property
    def state(self) -> NormalState:
        return self._normal_state


class NormalCoreModel(CoreModel):
    pass


class NormalProfile:
    def __init__(self, bg: NormalNullModel):
        self._multiple_hits: bool = True
        self._bg = bg

        alphabet = bg.state.alphabet
        emission_table = bg.state.emission_table()
        special_node = SpecialNode(
            S=MuteState("S", alphabet),
            N=NormalState("N", alphabet, emission_table),
            B=MuteState("B", alphabet),
            E=MuteState("E", alphabet),
            J=NormalState("J", alphabet, emission_table),
            C=NormalState("C", alphabet, emission_table),
            T=MuteState("T", alphabet),
        )
        hmm = HMM(alphabet)
        hmm.add_state(special_node.S, LOG1)
        hmm.add_state(special_node.N)
        hmm.add_state(special_node.B)
        hmm.add_state(special_node.E)
        hmm.add_state(special_node.J)
        hmm.add_state(special_node.C)
        hmm.add_state(special_node.T)

        self._hmm = hmm
        self._special_node = special_node
        self._special_trans = SpecialTrans()

        self._core_nodes: List[Node] = []

    def core_model(self):
        return NormalCoreModel(self._hmm, self._core_nodes, self._finalize)

    @property
    def length(self):
        return len(self._core_nodes)

    def set_multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def lr(self, seq: bytes) -> HMMERResult:
        self._set_target_length(seq)
        score0 = self._bg.likelihood(seq)
        result = self._viterbi(seq)
        score = result.score - score0
        return HMMERResult(score, seq, result.path)

    def _finalize(self):
        self._set_fragment_length()

    def _set_fragment_length(self):
        if self.length == 0:
            return

        B = self._special_node.B
        E = self._special_node.E

        # Uniform local alignment fragment length distribution
        t = self._special_trans
        t.BM = log(2) - log(self.length) - log(self.length + 1)
        t.ME = 0.0
        for node in self._core_nodes:
            self._hmm.set_trans(B, node.M, t.BM)
            self._hmm.set_trans(node.M, E, t.ME)

        for node in self._core_nodes[1:]:
            self._hmm.set_trans(node.D, E, 0.0)

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

        t = self._special_trans

        t.NN = t.CC = t.JJ = lp
        t.NB = t.CT = t.JB = l1p
        t.RR = lr
        t.EC = t.EJ = lq

        node = self._special_node

        self._hmm.set_trans(node.S, node.B, t.NB)
        self._hmm.set_trans(node.S, node.N, t.NN)
        self._hmm.set_trans(node.N, node.N, t.NN)
        self._hmm.set_trans(node.N, node.B, t.NB)

        self._hmm.set_trans(node.E, node.T, t.EC + t.CT)
        self._hmm.set_trans(node.E, node.C, t.EC + t.CC)
        self._hmm.set_trans(node.C, node.C, t.CC)
        self._hmm.set_trans(node.C, node.T, t.CT)

        self._hmm.set_trans(node.E, node.B, t.EJ + t.JB)
        self._hmm.set_trans(node.E, node.J, t.EJ + t.JJ)
        self._hmm.set_trans(node.J, node.J, t.JJ)
        self._hmm.set_trans(node.J, node.B, t.JB)

        self._bg.set_trans(t.RR)

    def _viterbi(self, seq: bytes) -> PathScore:
        self._set_target_length(seq)
        return self._hmm.viterbi(seq, self._special_node.T)


def create_hmmer_profile(reader: HMMEReader) -> NormalProfile:

    alphabet = Alphabet(reader.alphabet)
    R = NormalState("R", alphabet, reader.insert(0))
    R.normalize()
    hmmer = NormalProfile(NormalNullModel(R))

    with hmmer.core_model() as core:
        for m in range(1, reader.M + 1):
            node = Node(
                M=NormalState(f"M{m}", alphabet, reader.match(m)),
                I=NormalState(f"I{m}", alphabet, reader.insert(m)),
                D=MuteState(f"D{m}", alphabet),
            )
            node.M.normalize()
            node.I.normalize()
            trans = Trans(**reader.trans(m - 1))
            trans.normalize()
            core.add_node(node, trans)

    return hmmer
