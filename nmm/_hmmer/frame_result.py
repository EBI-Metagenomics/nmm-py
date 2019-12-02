from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval, SearchResult
from .frame_core import FramePath, FrameStep
from .._log import LOG1
from .._state import MuteState, FrameState, CodonState
from .._codon import Codon
from .codon import CodonFragment, CodonPath


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

    def decode_codons(self) -> CodonFragment:
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
    def interval(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def score(self) -> float:
        return self._score


def _create_path(steps: List[FrameStep]):
    path = FramePath()
    for step in steps:
        path.append_frame_step(step.state, step.seq_len)
    return path
