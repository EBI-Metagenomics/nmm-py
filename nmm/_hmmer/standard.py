import pathlib
from io import TextIOBase
from math import log
from typing import List, NamedTuple, Union

import hmmer_reader

from .._alphabet import Alphabet
from .._hmm import HMM
from .._log import LOG0
from .._state import MuteState, NormalState
from .background import BackgroundModel
from .core import HMMERCoreModel
from .._hmm import PathScore
from .path import HMMERResult
from .core import SpecialTrans, Trans


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


class HMMERProfile:
    def __init__(self, bg_model: BackgroundModel):
        self._multiple_hits: bool = True
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

        self._core_nodes: List[Node] = []

    def core_model(self):
        return HMMERCoreModel(self._hmm, self._core_nodes, self._finalize)

    @property
    def length(self):
        return len(self._core_nodes)

    def set_multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def lr(self, seq: str) -> HMMERResult:
        self._set_target_length(seq)
        score0 = self._bg_model.likelihood(seq)
        result = self._viterbi(seq)
        score = result.score - score0
        return HMMERResult(score, seq.encode(), result.path)

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

    def _set_target_length(self, seq: str):
        from math import exp

        L = len(seq.encode())
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

        self._bg_model.set_trans(t.RR)

    def _viterbi(self, seq: str) -> PathScore:
        self._set_target_length(seq)
        return self._hmm.viterbi(seq, self._special_node.T)


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
