import pathlib
from io import TextIOBase
from math import log
from typing import List, NamedTuple, Union, Dict
from dataclasses import dataclass

import hmmer_reader

from ._alphabet import Alphabet
from ._hmm import HMM
from ._log import LOG0
from ._state import MuteState, NormalState
from ._path import Path


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


@dataclass
class Trans:
    MM: float = LOG0
    MI: float = LOG0
    MD: float = LOG0
    IM: float = LOG0
    II: float = LOG0
    DM: float = LOG0
    DD: float = LOG0

    def normalize(self):
        from numpy import logaddexp

        m_norm: float = logaddexp(logaddexp(self.MM, self.MI), self.MD)
        self.MM -= m_norm
        self.MI -= m_norm
        self.MD -= m_norm

        i_norm: float = logaddexp(self.IM, self.II)
        self.IM -= i_norm
        self.II -= i_norm

        d_norm: float = logaddexp(self.DM, self.DD)
        self.DM -= d_norm
        self.DD -= d_norm


@dataclass
class SpecialTrans:
    NN: float = 0.0
    NB: float = 0.0
    EC: float = 0.0
    CC: float = 0.0
    CT: float = 0.0
    EJ: float = 0.0
    JJ: float = 0.0
    JB: float = 0.0


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


class HMMERProfile:
    def __init__(self, bg_model: BackgroundModel):
        self._alphabet = bg_model.state.alphabet
        self._bg_model = bg_model
        emission_table = bg_model.state.emission_table()
        self._hmm = HMM(self._alphabet)

        self._special_node = SpecialNode(
            S=MuteState("S", self._alphabet),
            N=NormalState("N", self._alphabet, emission_table),
            B=MuteState("B", self._alphabet),
            E=MuteState("E", self._alphabet),
            J=NormalState("J", self._alphabet, emission_table),
            C=NormalState("C", self._alphabet, emission_table),
            T=MuteState("T", self._alphabet),
        )
        self._hmm.add_state(self._special_node.S, 0.0)
        self._hmm.add_state(self._special_node.N)
        self._hmm.add_state(self._special_node.B)
        self._hmm.add_state(self._special_node.E)
        self._hmm.add_state(self._special_node.J)
        self._hmm.add_state(self._special_node.C)
        self._hmm.add_state(self._special_node.T)

        self._special_trans = SpecialTrans()
        self._assign_special_trans()

        self._core_nodes: List[Node] = []

    def _assign_special_trans(self):
        node = self._special_node
        trans = self._special_trans

        self._hmm.set_trans(node.S, node.B, trans.NB)
        self._hmm.set_trans(node.S, node.N, trans.NN)
        self._hmm.set_trans(node.N, node.N, trans.NN)
        self._hmm.set_trans(node.N, node.B, trans.NB)

        self._hmm.set_trans(node.E, node.T, trans.EC + trans.CT)
        self._hmm.set_trans(node.E, node.C, trans.EC + trans.CC)
        self._hmm.set_trans(node.C, node.C, trans.CC)
        self._hmm.set_trans(node.C, node.T, trans.CT)

        self._hmm.set_trans(node.E, node.B, trans.EJ + trans.JB)
        self._hmm.set_trans(node.E, node.J, trans.EJ + trans.JJ)
        self._hmm.set_trans(node.J, node.J, trans.JJ)
        self._hmm.set_trans(node.J, node.B, trans.JB)

    def core_model(self):
        return HMMERCoreModel(self._hmm, self._core_nodes, self._finalize)

    @property
    def length(self):
        return len(self._core_nodes)

    def _finalize(self):
        if self.length == 0:
            return
        B = self._special_node.B
        E = self._special_node.E
        # Uniform local alignment fragment length distribution
        t_BM = log(2) - log(self.length) - log(self.length + 1)
        t_ME = 0.0
        for node in self._core_nodes:
            self._hmm.set_trans(B, node.M, t_BM)
            self._hmm.set_trans(node.M, E, t_ME)

        for node in self._core_nodes[1:]:
            self._hmm.set_trans(node.D, E, 0.0)

    def _set_target_length(self, seq: str):
        L = len(seq.encode())
        p = log(L) - log(L + 2)
        q = LOG0
        r = log(L) - log(L + 1)

        special = self._special_states
        self._hmm.set_trans(self._special_node.S, self._special_node.N, p)
        self._hmm.set_trans(self._special_node.N, self._special_node.N, p)

        self._hmm.set_trans(self._special_node.E, self._special_node.C, p)
        self._hmm.set_trans(self._special_node.C, self._special_node.C, p)

        self._hmm.set_trans(self._special_node.J, self._special_node.J, p)

        self._hmm.set_trans(self._special_node.E, self._special_node.J, p * q)
        self._hmm.set_trans(self._special_node.E, self._special_node.B, q)
        self._bg_model.set_trans(r)

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def viterbi(self, seq: str) -> float:
        self._set_target_length(seq)
        return self._hmm.viterbi(seq, self._special_states["T"])

    def lr(self, seq: str) -> float:
        self._set_target_length(seq)
        score0 = self._bg_model.likelihood(seq)
        score1 = self.viterbi(seq)
        return score1 - score0


class HMMERCoreModel:
    def __init__(self, hmm: HMM, core_nodes: List[Node], finalize):
        self._hmm = hmm
        self._core_nodes = core_nodes
        self._finalize = finalize

    def add_node(self, node: Node, trans: Trans):
        self._hmm.add_state(node.M)
        self._hmm.add_state(node.I)
        self._hmm.add_state(node.D)

        self._core_nodes.append(node)

        if len(self._core_nodes) == 1:
            return

        prev = self._core_nodes[-2]
        self._hmm.set_trans(prev.M, node.M, trans.MM)
        self._hmm.set_trans(prev.M, prev.I, trans.MI)
        self._hmm.set_trans(prev.M, node.D, trans.MD)
        self._hmm.set_trans(prev.I, node.M, trans.IM)
        self._hmm.set_trans(prev.I, prev.I, trans.II)
        self._hmm.set_trans(prev.D, node.M, trans.DM)
        self._hmm.set_trans(prev.D, node.D, trans.DD)

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


def read_hmmer(file: Union[str, pathlib.Path, TextIOBase]) -> HMMERProfile:
    if isinstance(file, str):
        file = pathlib.Path(file)

    if isinstance(file, pathlib.Path):
        if not file.exists():
            raise ValueError(f"`{file}` does not exist.")

        if not file.is_file():
            raise ValueError(f"`{file}` is not a file.")

    hmmfile = hmmer_reader.read(file)
    alphabet = Alphabet(hmmfile.alphabet)
    R = NormalState("R", alphabet, hmmfile.insert(0))
    R.normalize()
    hmmer = HMMERProfile(BackgroundModel(R))

    with hmmer.core_model() as core:
        for m in range(1, hmmfile.M + 1):
            node = Node(
                M=NormalState(f"M{m}", alphabet, hmmfile.match(m)),
                I=NormalState(f"I{m}", alphabet, hmmfile.insert(m)),
                D=MuteState(f"D{m}", alphabet),
            )
            node.M.normalize()
            node.I.normalize()
            trans = Trans(**hmmfile.trans(m - 1))
            trans.normalize()
            core.add_node(node, trans)

    return hmmer
