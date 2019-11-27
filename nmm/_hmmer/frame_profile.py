from math import log
from typing import Any, Dict, List, NamedTuple, Tuple

from hmmer_reader import HMMEReader

from .._alphabet import Alphabet
from .._base import Base
from .._codon import CodonTable
from .._gencode import GeneticCode
from .._hmm import HMM
from .._log import LOG0, LOG1
from .._path import CPath, Path
from .._state import FrameState, MuteState
from .core import CoreModel, NullModel
from .result import Result
from .transition import SpecialTransitions, Transitions

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


class FrameNullModel(NullModel):
    def __init__(self, state: FrameState):
        super().__init__(state)
        self._frame_state = state

    @property
    def state(self) -> FrameState:
        return self._frame_state


class FrameCoreModel(CoreModel):
    pass


class FrameProfile:
    def __init__(self, bg: FrameNullModel):
        self._multiple_hits: bool = True
        self._bg = bg

        alphabet = bg.state.alphabet
        base = bg.state.base
        codon = bg.state.codon
        epsilon = bg.state.epsilon
        special_node = SpecialNode(
            S=MuteState(b"S", alphabet),
            N=FrameState(b"N", base, codon, epsilon),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=FrameState(b"J", base, codon, epsilon),
            C=FrameState(b"C", base, codon, epsilon),
            T=MuteState(b"T", alphabet),
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
        self._special_trans = SpecialTransitions()

        self._core_nodes: List[Node] = []

    def core_model(self):
        return FrameCoreModel(self._hmm, self._core_nodes, self._finalize)

    @property
    def length(self):
        return len(self._core_nodes)

    def set_multiple_hits(self, multiple_hits: bool):
        self._multiple_hits = multiple_hits

    @property
    def hmm(self) -> HMM:
        return self._hmm

    def lr(self, seq: bytes) -> Tuple[Result, Result]:
        self._set_target_length(seq)
        score0 = self._bg.likelihood(seq)
        score1, path = self._viterbi(seq)
        score = score1 - score0
        codon_seq, codon_path = self._convert_to_codon_path(seq, path)
        return (
            Result(score, seq, path),
            Result(score, codon_seq, codon_path),
        )

    def _convert_to_codon_path(self, seq: bytes, path):
        nseq: List[bytes] = []
        npath = Path()
        start: int = 0
        for step in path.steps():
            state = self._hmm.states()[step.state.imm_state]
            if step.seq_len == 0:
                npath.append(state, 0)
            else:
                fstate: FrameState = state
                decoded_codon = fstate.decode(seq[start : start + step.seq_len])
                nseq.append(decoded_codon.codon)
                npath.append(fstate, 3)
            start += step.seq_len

        return (b"".join(nseq), npath)

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
            self._hmm.set_transition(B, node.M, t.BM)
            self._hmm.set_transition(node.M, E, t.ME)

        for node in self._core_nodes[1:]:
            self._hmm.set_transition(node.D, E, 0.0)

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

        self._hmm.set_transition(node.S, node.B, t.NB)
        self._hmm.set_transition(node.S, node.N, t.NN)
        self._hmm.set_transition(node.N, node.N, t.NN)
        self._hmm.set_transition(node.N, node.B, t.NB)

        self._hmm.set_transition(node.E, node.T, t.EC + t.CT)
        self._hmm.set_transition(node.E, node.C, t.EC + t.CC)
        self._hmm.set_transition(node.C, node.C, t.CC)
        self._hmm.set_transition(node.C, node.T, t.CT)

        self._hmm.set_transition(node.E, node.B, t.EJ + t.JB)
        self._hmm.set_transition(node.E, node.J, t.EJ + t.JJ)
        self._hmm.set_transition(node.J, node.J, t.JJ)
        self._hmm.set_transition(node.J, node.B, t.JB)

        self._bg.set_transition(t.RR)

    def _viterbi(self, seq: bytes) -> Tuple[float, CPath]:
        self._set_target_length(seq)
        return self._hmm.viterbi(seq, self._special_node.T)


def _infer_codon_lprobs(aa_lprobs: Dict[bytes, float], gencode: GeneticCode):
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

    lprobs: Dict[bytes, list] = {bytes([base]): [] for base in alphabet.symbols}
    lprob_norm = log(3)
    for codon, lprob in codon_lprobs.items():
        lprobs[codon[0:1]] += [lprob - lprob_norm]
        lprobs[codon[1:2]] += [lprob - lprob_norm]
        lprobs[codon[2:3]] += [lprob - lprob_norm]

    return {b: logsumexp(lp) for b, lp in lprobs.items()}


class _FrameStateFactory:
    def __init__(self, bases: Alphabet, gcode: GeneticCode, epsilon: float):
        self._bases = bases
        self._gcode = gcode
        self._epsilon = epsilon

    def create(self, name: bytes, aa_lprobs: Dict[bytes, float]) -> FrameState:
        codon_lprobs = _infer_codon_lprobs(aa_lprobs, self._gcode)
        base_lprobs = _infer_base_lprobs(codon_lprobs, self._bases)
        base = Base(self._bases, base_lprobs)
        codon = CodonTable(self._bases, codon_lprobs)
        return FrameState(name, base, codon, self._epsilon)


def create_frame_profile(reader: HMMEReader, epsilon: float = 0.1) -> FrameProfile:
    bases = Alphabet(b"ACGU")
    ffact = _FrameStateFactory(bases, GeneticCode(), epsilon)
    R = ffact.create(b"R", _bytes_dict(reader.insert(0)))

    # TODO: the null model is not property set.
    # It is supposed to be temporary.
    hmmer = FrameProfile(FrameNullModel(R))

    with hmmer.core_model() as core:
        for m in range(1, reader.M + 1):
            node = Node(
                M=ffact.create(f"M{m}".encode(), _bytes_dict(reader.match(m))),
                I=ffact.create(f"I{m}".encode(), _bytes_dict(reader.insert(m))),
                D=MuteState(f"D{m}".encode(), bases),
            )
            trans = Transitions(**reader.trans(m - 1))
            trans.normalize()
            core.add_node(node, trans)

    return hmmer


def _bytes_dict(d: Dict[str, Any]):
    return {k.encode(): v for k, v in d.items()}
