from typing import Iterator, List, Sequence, Tuple, TypeVar, Union

from .._codon import Codon
from .._ffi import ffi
from .._log import LOG1
from .._path import CPath
from .._state import CodonState, FrameState, MuteState
from .._step import CStep
from .codon import CodonFragment, CodonPath, CodonSearchResult
from .result import Fragment, Interval, SearchResult


class FrameStep(CStep):
    def __init__(self, imm_step: ffi.CData, state: Union[MuteState, FrameState]):
        super().__init__(imm_step)
        self._state = state

    @property
    def state(self) -> Union[MuteState, FrameState]:
        return self._state


T = TypeVar("T", bound="FramePath")


class FramePath(CPath):
    def __init__(self):
        super().__init__()
        self._steps: List[FrameStep] = []

    def append_frame_step(
        self, state: Union[MuteState, FrameState], seq_len: int
    ) -> FrameStep:
        cstep = self.append_cstep(state, seq_len)
        self._steps.append(FrameStep(cstep.imm_step, state))
        return self._steps[-1]

    def steps(self) -> Iterator[FrameStep]:
        return iter(self._steps)


class FrameFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: FramePath, homologous: bool,
    ):
        super().__init__(homologous)
        self._path = path
        self._sequence = sequence

    @property
    def sequence(self) -> bytes:
        return self._sequence

    def items(self) -> Iterator[Tuple[bytes, FrameStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def decode(self) -> CodonFragment:
        nseq: List[Codon] = []
        npath = CodonPath()

        start: int = 0
        seq = self.sequence
        for step in self._path.steps():
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, step.state.alphabet)
                npath.append_codon_step(mstate, 0)
            else:
                assert isinstance(step.state, FrameState)

                codon = step.state.decode(seq[start : start + step.seq_len])[0]
                nseq.append(codon)

                cstate = CodonState(step.state.name, step.state.alphabet, {codon: LOG1})
                npath.append_codon_step(cstate, 3)

            start += step.seq_len

        return CodonFragment(nseq, npath, self.homologous)

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class FrameSearchResult(SearchResult):
    def __init__(self, score: float, sequence: bytes, path: FramePath):
        self._score = score

        self._fragments: List[FrameFragment] = []
        self._intervals: List[Interval] = []

        steps = list(path.steps())
        for fragi, stepi, homologous in self._create_fragments(path):
            spath = _create_path(steps[stepi.start : stepi.end])
            seq = sequence[fragi.start : fragi.end]
            frag = FrameFragment(seq, spath, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> Sequence[FrameFragment]:
        return self._fragments

    @property
    def intervals(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def score(self) -> float:
        return self._score

    def decode(self) -> CodonSearchResult:
        fragments: List[CodonFragment] = []
        intervals: List[Interval] = []

        start = end = 0
        for i, frag in enumerate(self._fragments):

            codon_frag = frag.decode()
            end += len(codon_frag.sequence)

            fragments.append(codon_frag)
            intervals.append(Interval(start, end))

            start = end

        return CodonSearchResult(self.score, fragments, intervals)


def _create_path(steps: List[FrameStep]):
    path = FramePath()
    for step in steps:
        path.append_frame_step(step.state, step.seq_len)
    return path
