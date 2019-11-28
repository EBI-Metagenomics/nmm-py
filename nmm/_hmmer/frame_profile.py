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
        self, bases: Alphabet, gcode: GeneticCode, epsilon: float,
    ):
        self._bases = bases
        self._gcode = gcode
        self._epsilon = epsilon

    def create(self, name: bytes, aa_lprobs: Dict[bytes, float]) -> FrameState:
        codon_lprobs = _infer_codon_lprobs(aa_lprobs, self._gcode)
        base_lprobs = _infer_base_lprobs(codon_lprobs, self._bases)
        base_table = Base(self._bases, base_lprobs)
        codon_table = CodonTable(self._bases, codon_lprobs)
        return FrameState(name, base_table, codon_table, self._epsilon)

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
        aa_lprobs: Dict[bytes, float],
        nodes_trans: Sequence[Tuple[FrameNode, Transitions]],
    ):
        super().__init__()

        R = fstate_factory.create(b"R", aa_lprobs)
        self._null_model = FrameNullModel(R)

        special_node = FrameSpecialNode(
            S=MuteState(b"S", fstate_factory.bases),
            N=fstate_factory.create(b"N", aa_lprobs),
            B=MuteState(b"B", fstate_factory.bases),
            E=MuteState(b"E", fstate_factory.bases),
            J=fstate_factory.create(b"J", aa_lprobs),
            C=fstate_factory.create(b"C", aa_lprobs),
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


def create_frame_profile(reader: HMMEReader, epsilon: float = 0.1) -> FrameProfile:

    bases = Alphabet(b"ACGU")
    null_lprobs = _dict(reader.insert(0))
    ffact = FrameStateFactory(bases, GeneticCode(), epsilon)

    nodes_trans: List[Tuple[FrameNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        M = ffact.create(f"M{m}".encode(), _dict(reader.match(m)))
        I = ffact.create(f"I{m}".encode(), _dict(reader.insert(m)))
        D = MuteState(f"D{m}".encode(), bases)

        node = FrameNode(M, I, D,)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return FrameProfile(ffact, null_lprobs, nodes_trans)


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


def _dict(d: Dict[str, Any]):
    return {k.encode(): v for k, v in d.items()}
