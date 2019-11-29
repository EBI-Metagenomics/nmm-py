from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval, SearchResult
from .frame_core import FramePath, FrameStep
from .._state import MuteState, FrameState


class FrameFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: FramePath, homologous: bool,
    ):
        super().__init__(sequence, homologous)
        self._path = path

    def items(self) -> Iterator[Tuple[bytes, FrameStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def convert_to_codons(self):
        nseq: List[bytes] = []
        npath = FramePath()

        start: int = 0
        seq = self.sequence
        for step in self._path.steps():
            if isinstance(step.state, MuteState):
                npath.append_frame_step(step.state, 0)
            else:
                assert isinstance(step.state, FrameState)
                decoded_codon = step.state.decode(seq[start : start + step.seq_len])
                nseq.append(decoded_codon.codon)
                npath.append_frame_step(step.state, 3)

            start += step.seq_len

        return FrameFragment(b"".join(nseq), npath, self.homologous)

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
