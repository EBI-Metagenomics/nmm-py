from typing import Any, Dict, List, Sequence

from hmmer_reader import HMMEReader

from .._alphabet import Alphabet
from .._state import MuteState, NormalState
from .core import Profile
from .result import Result
from .transition import Transitions
from .standard_core import (
    StandardNullModel,
    StandardCoreModel,
    StandardSpecialNode,
    StandardNode,
)


class StandardProfile(Profile):
    def __init__(self, null_model: StandardNullModel):
        self._null_model = null_model

        alphabet = null_model.state.alphabet
        emission_table = null_model.state.emission_table()
        special_node = StandardSpecialNode(
            S=MuteState(b"S", alphabet),
            N=NormalState(b"N", alphabet, emission_table),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=NormalState(b"J", alphabet, emission_table),
            C=NormalState(b"C", alphabet, emission_table),
            T=MuteState(b"T", alphabet),
        )
        super().__init__(special_node)

        self._special_node = special_node
        self._core_nodes: List[StandardNode] = []

    @property
    def null_model(self) -> StandardNullModel:
        return self._null_model

    @property
    def special_node(self) -> StandardSpecialNode:
        return self._special_node

    @property
    def length(self):
        return len(self._core_nodes)

    def core_nodes(self) -> Sequence[StandardNode]:
        return self._core_nodes

    @property
    def core_model(self):
        return StandardCoreModel(self._hmm, self._core_nodes, self._finalize)

    def lr(self, seq: bytes) -> Result:
        self._set_target_length(seq)
        score0 = self.null_model.likelihood(seq)
        score1, path = self._viterbi(seq)
        score = score1 - score0
        return Result(score, seq, path)


def create_hmmer_profile(reader: HMMEReader) -> StandardProfile:

    alphabet = Alphabet(reader.alphabet.encode())
    # TODO: the null model is not property set.
    # It is supposed to be temporary.
    R = NormalState(b"R", alphabet, _dict(reader.insert(0)))
    R.normalize()
    hmmer = StandardProfile(StandardNullModel(R))

    with hmmer.core_model as core:
        for m in range(1, reader.M + 1):
            node = StandardNode(
                M=NormalState(f"M{m}".encode(), alphabet, _dict(reader.match(m))),
                I=NormalState(f"I{m}".encode(), alphabet, _dict(reader.insert(m))),
                D=MuteState(f"D{m}".encode(), alphabet),
            )
            node.M.normalize()
            node.I.normalize()
            trans = Transitions(**reader.trans(m - 1))
            trans.normalize()
            core.add_node(node, trans)

    return hmmer


def _dict(d: Dict[str, Any]):
    return {k.encode(): v for k, v in d.items()}
