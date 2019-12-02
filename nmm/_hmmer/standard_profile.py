from typing import Any, Dict, List, Sequence, Tuple

from hmmer_reader import HMMEReader

from .._alphabet import Alphabet
from .._state import MuteState, NormalState
from .model import Transitions
from .profile import Profile
from .standard import StandardSearchResult
from .standard_model import (
    StandardAltModel,
    StandardNode,
    StandardNullModel,
    StandardSpecialNode,
)


class StandardProfile(Profile):
    def __init__(
        self,
        alphabet: Alphabet,
        null_lprobs: Dict[bytes, float],
        nodes_trans: Sequence[Tuple[StandardNode, Transitions]],
    ):
        super().__init__()

        R = NormalState(b"R", alphabet, null_lprobs)
        R.normalize()
        self._null_model = StandardNullModel(R)

        emission_table = R.emission_table()
        special_node = StandardSpecialNode(
            S=MuteState(b"S", alphabet),
            N=NormalState(b"N", alphabet, emission_table),
            B=MuteState(b"B", alphabet),
            E=MuteState(b"E", alphabet),
            J=NormalState(b"J", alphabet, emission_table),
            C=NormalState(b"C", alphabet, emission_table),
            T=MuteState(b"T", alphabet),
        )

        self._alt_model = StandardAltModel(special_node, nodes_trans)
        self._set_fragment_length()

    @property
    def null_model(self) -> StandardNullModel:
        return self._null_model

    @property
    def alt_model(self) -> StandardAltModel:
        return self._alt_model

    def search(self, seq: bytes) -> StandardSearchResult:
        self._set_target_length(len(seq))
        score0 = self.null_model.likelihood(seq)
        score1, path = self.alt_model.viterbi(seq)
        score = score1 - score0
        return StandardSearchResult(score, seq, path)


def create_standard_profile(reader: HMMEReader) -> StandardProfile:

    alphabet = Alphabet(reader.alphabet.encode())
    null_lprobs = _dict(reader.insert(0))

    nodes_trans: List[Tuple[StandardNode, Transitions]] = []

    for m in range(1, reader.M + 1):
        M = NormalState(f"M{m}".encode(), alphabet, _dict(reader.match(m)))
        M.normalize()

        I = NormalState(f"I{m}".encode(), alphabet, _dict(reader.insert(m)))
        I.normalize()

        D = MuteState(f"D{m}".encode(), alphabet)

        node = StandardNode(M, I, D)

        trans = Transitions(**reader.trans(m - 1))
        trans.normalize()

        nodes_trans.append((node, trans))

    return StandardProfile(alphabet, null_lprobs, nodes_trans)


def _dict(d: Dict[str, Any]):
    return {k.encode(): v for k, v in d.items()}
