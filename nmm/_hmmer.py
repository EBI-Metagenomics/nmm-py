import pathlib
from io import TextIOBase
from math import log
from typing import List, NamedTuple, Union, Dict

import hmmer_reader

from ._alphabet import Alphabet
from ._hmm import HMM
from ._log import LOG0
from ._state import MuteState, NormalState
from ._path import Path


_Node = NamedTuple("Node", [("M", NormalState), ("I", NormalState), ("D", MuteState)])


class _BackgroundModel:
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
    def __init__(self, bg_model: _BackgroundModel):
        self._alphabet = bg_model.state.alphabet
        self._bg_model = bg_model
        emission_table = bg_model.state.emission_table()

        self._special_states = {
            "S": MuteState("S", self._alphabet),
            "N": NormalState("N", self._alphabet, emission_table),
            "B": MuteState("B", self._alphabet),
            "E": MuteState("E", self._alphabet),
            "J": NormalState("J", self._alphabet, emission_table),
            "C": NormalState("C", self._alphabet, emission_table),
            "T": MuteState("T", self._alphabet),
        }
        self._hmm = HMM(self._alphabet)

        special = self._special_states
        self._hmm.add_state(special["S"], 0.0)
        self._hmm.add_state(special["N"])
        self._hmm.add_state(special["B"])
        self._hmm.add_state(special["E"])
        self._hmm.add_state(special["J"])
        self._hmm.add_state(special["C"])
        self._hmm.add_state(special["T"])

        self._hmm.set_trans(special["S"], special["B"], 0.0)
        self._hmm.set_trans(special["N"], special["B"], 0.0)

        self._hmm.set_trans(special["E"], special["T"], 0.0)
        self._hmm.set_trans(special["C"], special["T"], 0.0)

        self._hmm.set_trans(special["J"], special["B"], 0.0)

        self._core_nodes: List[_Node] = []

    def core_model(self):
        return _HMMERCoreModel(self, self._core_nodes, self._finalize)

    def _finalize(self):
        B = self._special_states["B"]
        E = self._special_states["E"]
        for node in self._core_nodes:
            self._hmm.set_trans(B, node.M, 0.0)
            self._hmm.set_trans(node.M, E, 0.0)

        for node in self._core_nodes[1:]:
            self._hmm.set_trans(node.D, E, 0.0)

    def _set_target_length(self, seq: str):
        L = len(seq.encode())
        p = log(L) - log(L + 2)
        q = LOG0
        r = log(L) - log(L + 1)

        special = self._special_states
        self._hmm.set_trans(special["S"], special["N"], p)
        self._hmm.set_trans(special["N"], special["N"], p)

        self._hmm.set_trans(special["E"], special["C"], p)
        self._hmm.set_trans(special["C"], special["C"], p)

        self._hmm.set_trans(special["J"], special["J"], p)

        self._hmm.set_trans(special["E"], special["J"], p * q)
        self._hmm.set_trans(special["E"], special["B"], q)
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


class _HMMERCoreModel:
    def __init__(self, hmmer: HMMERProfile, core_nodes: List[_Node], finalize):
        self._hmmer = hmmer
        self._core_nodes = core_nodes
        self._finalize = finalize

    def add_node(
        self, M: NormalState, I: NormalState, D: MuteState, trans: Dict[str, float]
    ):
        self._hmmer.hmm.add_state(M)
        self._hmmer.hmm.add_state(I)
        self._hmmer.hmm.add_state(D)

        self._core_nodes.append(_Node(M=M, I=I, D=D))

        if len(self._core_nodes) == 1:
            return

        prev = self._core_nodes[-2]
        hmm = self._hmmer.hmm
        hmm.set_trans(prev.M, M, trans["MM"])
        hmm.set_trans(prev.M, prev.I, trans["MI"])
        hmm.set_trans(prev.M, D, trans["MD"])
        hmm.set_trans(prev.I, M, trans["IM"])
        hmm.set_trans(prev.I, prev.I, trans["II"])
        hmm.set_trans(prev.D, M, trans["DM"])
        hmm.set_trans(prev.D, D, trans["DD"])
        hmm.normalize_trans(prev.M)
        hmm.normalize_trans(prev.I)
        hmm.normalize_trans(prev.D)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if len(self._core_nodes) > 0:
            self._hmmer.hmm.del_state(self._core_nodes[0].D)
            self._hmmer.hmm.del_state(self._core_nodes[-1].I)
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
    R = NormalState("R", alphabet, hmmfile.insert(0, True))
    R.normalize()
    bg = _BackgroundModel(R)
    hmmer = HMMERProfile(bg)

    with hmmer.core_model() as core:
        for m in range(1, hmmfile.M + 1):
            M = NormalState(f"M{m}", alphabet, hmmfile.match(m, True))
            M.normalize()
            I = NormalState(f"I{m}", alphabet, hmmfile.insert(m, True))
            I.normalize()
            D = MuteState(f"D{m}", alphabet)
            trans = hmmfile.trans(m - 1, True)
            core.add_node(M, I, D, trans)

    return hmmer
