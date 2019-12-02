from typing import Iterator, List, Tuple, TypeVar, Union, Sequence

from .._ffi import ffi
from .._path import CPath
from .._state import MuteState, NormalState
from .._step import CStep
from .result import Fragment, SearchResult, Interval


class AminoAcidStep(CStep):
    def __init__(self, imm_step: ffi.CData, state: Union[MuteState, NormalState]):
        super().__init__(imm_step)
        self._state = state

    @property
    def state(self) -> Union[MuteState, NormalState]:
        return self._state


T = TypeVar("T", bound="AminoAcidPath")


class AminoAcidPath(CPath):
    def __init__(self):
        super().__init__()
        self._steps: List[AminoAcidStep] = []

    def append_amino_acid_step(
        self, state: Union[MuteState, NormalState], seq_len: int
    ) -> AminoAcidStep:
        cstep = self.append_cstep(state, seq_len)
        self._steps.append(AminoAcidStep(cstep.imm_step, state))
        return self._steps[-1]

    def steps(self) -> Iterator[AminoAcidStep]:
        return iter(self._steps)


class AminoAcidFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: AminoAcidPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> bytes:
        return self._sequence

    def items(self) -> Iterator[Tuple[bytes, AminoAcidStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class AminoAcidSearchResult(SearchResult):
    def __init__(
        self,
        score: float,
        fragments: Sequence[AminoAcidFragment],
        intervals: Sequence[Interval],
    ):
        self._score = score

        self._fragments: List[AminoAcidFragment] = list(fragments)
        self._intervals: List[Interval] = list(intervals)

    @property
    def fragments(self) -> Sequence[AminoAcidFragment]:
        return self._fragments

    @property
    def intervals(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def score(self) -> float:
        return self._score
