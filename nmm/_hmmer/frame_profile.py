from math import log
from typing import Any, Dict, List, NamedTuple, Tuple, Sequence

from hmmer_reader import HMMEReader

from .._alphabet import Alphabet
from .._base import Base
from .._codon import CodonTable
from .._gencode import GeneticCode
from .._hmm import HMM
from .._log import LOG0, LOG1
from .._path import CPath, Path
from .._state import FrameState, MuteState
from .core import AltModel, NullModel
from .profile import Profile
from .result import SearchResult
from .frame_result import FrameSearchResult
from .transition import SpecialTransitions, Transitions
from .frame_core import (
    FrameAltModel,
    FrameNode,
    FrameNullModel,
    FrameSpecialNode,
)


class FrameStateFactory:
    def __init__(
        self,
        bases: Alphabet,
        gcode: GeneticCode,
        aa_lprobs: Dict[bytes, float],
        epsilon: float,
    ):
        self._bases = bases
        self._gcode = gcode
        self._epsilon = epsilon
        codon_lprobs = _infer_codon_lprobs(aa_lprobs, self._gcode)
        base_lprobs = _infer_base_lprobs(codon_lprobs, self._bases)
        self._base_table = Base(self._bases, base_lprobs)
        self._codon_table = CodonTable(self._bases, codon_lprobs)

    def create(self, name: bytes) -> FrameState:
        return FrameState(name, self._base_table, self._codon_table, self._epsilon)

    @property
    def bases(self) -> Alphabet:
        return self._bases

    @property
    def genetic_code(self) -> GeneticCode:
        return self._gcode

    @property
    def epsilon(self) -> float:
        return self._epsilon


class FrameProfile(Profile):
    def __init__(
        self,
        fstate_factory: FrameStateFactory,
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        super().__init__()

        R = fstate_factory.create(b"R")
        self._null_model = FrameNullModel(R)

        # emission_table = R.emission_table()
        # special_node = FrameSpecialNode(
        #     S=MuteState(b"S", alphabet),
        #     N=FrameState(b"N", alphabet, emission_table),
        #     B=MuteState(b"B", alphabet),
        #     E=MuteState(b"E", alphabet),
        #     J=FrameState(b"J", alphabet, emission_table),
        #     C=FrameState(b"C", alphabet, emission_table),
        #     T=MuteState(b"T", alphabet),
        # )

        # base = bg.state.base
        # codon = bg.state.codon
        # epsilon = bg.state.epsilon
        special_node = FrameSpecialNode(
            S=MuteState(b"S", fstate_factory.bases),
            N=fstate_factory.create(b"N"),
            B=MuteState(b"B", fstate_factory.bases),
            E=MuteState(b"E", fstate_factory.bases),
            J=fstate_factory.create(b"J"),
            C=fstate_factory.create(b"C"),
            T=MuteState(b"T", fstate_factory.bases),
        )

        self._alt_model = FrameAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> FrameNullModel:
        return self._null_model

    @property
    def alt_model(self) -> FrameAltModel:
        return self._alt_model

    def search(self, seq: bytes) -> FrameSearchResult:
        self._set_target_length(len(seq))
        score0 = self.null_model.likelihood(seq)
        score1, path = self.alt_model.viterbi(seq)
        score = score1 - score0
        return FrameSearchResult(score, seq, path)

    # def lr(self, seq: bytes) -> Tuple[SearchResult, SearchResult]:
    #     self._set_target_length(seq)
    #     score0 = self._bg.likelihood(seq)
    #     score1, path = self._viterbi(seq)
    #     score = score1 - score0
    #     codon_seq, codon_path = self._convert_to_codon_path(seq, path)
    #     return (
    #         SearchResult(score, seq, path),
    #         SearchResult(score, codon_seq, codon_path),
    #     )

    # def _convert_to_codon_path(self, seq: bytes, path):
    #     nseq: List[bytes] = []
    #     npath = Path()
    #     start: int = 0
    #     for step in path.steps():
    #         state = self._hmm.states()[step.state.imm_state]
    #         if step.seq_len == 0:
    #             npath.append(state, 0)
    #         else:
    #             fstate: FrameState = state
    #             decoded_codon = fstate.decode(seq[start : start + step.seq_len])
    #             nseq.append(decoded_codon.codon)
    #             npath.append(fstate, 3)
    #         start += step.seq_len

    #     return (b"".join(nseq), npath)

    # def _finalize(self):
    #     self._set_fragment_length()

    # def _set_fragment_length(self):
    #     if self.length == 0:
    #         return

    #     B = self._special_node.B
    #     E = self._special_node.E

    #     # Uniform local alignment fragment length distribution
    #     t = self._special_trans
    #     t.BM = log(2) - log(self.length) - log(self.length + 1)
    #     t.ME = 0.0
    #     for node in self._core_nodes:
    #         self._hmm.set_transition(B, node.M, t.BM)
    #         self._hmm.set_transition(node.M, E, t.ME)

    #     for node in self._core_nodes[1:]:
    #         self._hmm.set_transition(node.D, E, 0.0)

    # def _set_target_length(self, seq: bytes):
    #     from math import exp

    #     L = len(seq)
    #     if L == 0:
    #         return

    #     if self._multiple_hits:
    #         lq = -log(2)
    #     else:
    #         lq = LOG0

    #     q = exp(lq)
    #     lp = log(L) - log(L + 2 + q / (1 - q))
    #     l1p = log(2 + q / (1 - q)) - log(L + 2 + q / (1 - q))
    #     lr = log(L) - log(L + 1)

    #     t = self._special_trans

    #     t.NN = t.CC = t.JJ = lp
    #     t.NB = t.CT = t.JB = l1p
    #     t.RR = lr
    #     t.EC = t.EJ = lq

    #     node = self._special_node

    #     self._hmm.set_transition(node.S, node.B, t.NB)
    #     self._hmm.set_transition(node.S, node.N, t.NN)
    #     self._hmm.set_transition(node.N, node.N, t.NN)
    #     self._hmm.set_transition(node.N, node.B, t.NB)

    #     self._hmm.set_transition(node.E, node.T, t.EC + t.CT)
    #     self._hmm.set_transition(node.E, node.C, t.EC + t.CC)
    #     self._hmm.set_transition(node.C, node.C, t.CC)
    #     self._hmm.set_transition(node.C, node.T, t.CT)

    #     self._hmm.set_transition(node.E, node.B, t.EJ + t.JB)
    #     self._hmm.set_transition(node.E, node.J, t.EJ + t.JJ)
    #     self._hmm.set_transition(node.J, node.J, t.JJ)
    #     self._hmm.set_transition(node.J, node.B, t.JB)

    #     self._bg.set_transition(t.RR)


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


def create_frame_profile(reader: HMMEReader, epsilon: float = 0.1) -> FrameProfile:
    bases = Alphabet(b"ACGU")
    ffact = FrameStateFactory(bases, GeneticCode(), epsilon)
    R = ffact.create(b"R", _dict(reader.insert(0)))

    # TODO: the null model is not property set.
    # It is supposed to be temporary.
    hmmer = FrameProfile(FrameNullModel(R))

    with hmmer.core_model() as core:
        for m in range(1, reader.M + 1):
            node = Node(
                M=ffact.create(f"M{m}".encode(), _dict(reader.match(m))),
                I=ffact.create(f"I{m}".encode(), _dict(reader.insert(m))),
                D=MuteState(f"D{m}".encode(), bases),
            )
            trans = Transitions(**reader.trans(m - 1))
            trans.normalize()
            core.add_node(node, trans)

    return hmmer


def _dict(d: Dict[str, Any]):
    return {k.encode(): v for k, v in d.items()}
