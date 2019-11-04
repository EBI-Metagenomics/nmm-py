import pathlib
from io import TextIOBase
from math import log
from typing import List, NamedTuple, Union, Dict

from hmmer_reader import HMMEReader

from .._base import Base
from .._alphabet import Alphabet
from .._hmm import HMM, PathScore
from .._log import LOG0
from .._state import FrameState, MuteState
from .._gencode import GeneticCode
from .._codon import Codon
from .background import BackgroundModel
from .core import HMMERCoreModel, SpecialTrans, Trans
from .path import HMMERResult

Node = NamedTuple("Node", [("M", FrameState), ("I", FrameState), ("D", MuteState)])

SpecialNode = NamedTuple(
    "SpecialNode",
    [
        ("S", MuteState),
        ("N", FrameState),
        ("B", MuteState),
        ("E", MuteState),
        ("J", FrameState),
        ("C", FrameState),
        ("T", MuteState),
    ],
)


class HMMERFrameProfile:
    def __init__(self, bg_model: BackgroundModel):
        self._multiple_hits: bool = True
        self._alphabet = bg_model.state.alphabet
        self._bg_model = bg_model
        emission_table = bg_model.state.emission_table()
        self._hmm = HMM(self._alphabet)

        self._special_node = SpecialNode(
            S=MuteState("S", self._alphabet),
            N=FrameState("N", self._alphabet, emission_table),
            B=MuteState("B", self._alphabet),
            E=MuteState("E", self._alphabet),
            J=FrameState("J", self._alphabet, emission_table),
            C=FrameState("C", self._alphabet, emission_table),
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


def _infer_codon_lprobs(aa_lprobs: Dict[str, float], gencode: GeneticCode):
    from numpy import logaddexp

    codon_lprobs = []
    lprob_norm = LOG0
    for aa, lprob in aa_lprobs.items():

        codons = gencode.codons(aa)
        if len(codons) == 0:
            continue

        norm = log(len(codons))
        for codon in codons:
            codon_lprobs.append((codon, lprob - norm))
            lprob_norm = logaddexp(lprob_norm, codon_lprobs[-1][1])

    codon_lprobs = [(i[0], i[1] - lprob_norm) for i in codon_lprobs]
    return dict(codon_lprobs)


def _infer_base_lprobs(codon_lprobs, alphabet: Alphabet):
    from scipy.special import logsumexp

    lprobs: Dict[str, list] = {base: [] for base in alphabet.symbols}
    lprob_norm = log(3)
    for codon, lprob in codon_lprobs.items():
        lprobs[codon[0]] += [lprob - lprob_norm]
        lprobs[codon[1]] += [lprob - lprob_norm]
        lprobs[codon[2]] += [lprob - lprob_norm]

    return {b: logsumexp(lp) for b, lp in lprobs.items()}


def create_frame_profile(reader: HMMEReader) -> HMMERFrameProfile:
    gcode = GeneticCode()
    codon_lprobs = _infer_codon_lprobs(reader.compo, gcode)
    # codon_abc = Alphabet(reader.alphabet)

    bases_abc = Alphabet("ACGU")
    base_lprobs = _infer_base_lprobs(codon_lprobs, bases_abc)
    base = Base(bases_abc, base_lprobs)
    codon = Codon(bases_abc, codon_lprobs)

    epsilon = 0.1
    R = FrameState("R", base, codon, epsilon)
    # # M3 = FrameState("M3", base_emission, codon_emission, epsilon)
    # hmmer = HMMERProfile(BackgroundModel(R))

    # with hmmer.core_model() as core:
    #     for m in range(1, reader.M + 1):
    #         node = Node(
    #             M=NormalState(f"M{m}", alphabet, reader.match(m)),
    #             I=NormalState(f"I{m}", alphabet, reader.insert(m)),
    #             D=MuteState(f"D{m}", alphabet),
    #         )
    #         node.M.normalize()
    #         node.I.normalize()
    #         trans = Trans(**reader.trans(m - 1))
    #         trans.normalize()
    #         core.add_node(node, trans)

    # return hmmer
