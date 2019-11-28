from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval, SearchResult
from .frame_core import FramePath, FrameStep


class FrameFragment(Fragment):
    def __init__(
        self, seq: bytes, interval: Interval, path: FramePath, homologous: bool,
    ):
        super().__init__(seq, interval, homologous)
        self._path = path

    def items(self) -> Iterator[Tuple[bytes, FrameStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class FrameSearchResult(SearchResult):
    def __init__(self, score: float, seq: bytes, path: FramePath):
        self._score = score

        self._fragments: List[FrameFragment] = []

        steps = list(path.steps())
        for fragi, stepi, homologous in self._create_fragments(path):
            spath = _create_path(steps[stepi.start : stepi.end])
            self._fragments.append(FrameFragment(seq, fragi, spath, homologous))

    @property
    def fragments(self) -> Sequence[FrameFragment]:
        return self._fragments

    @property
    def score(self) -> float:
        return self._score


def _create_path(steps: List[FrameStep]):
    path = FramePath()
    for step in steps:
        path.append_frame_step(step.state, step.seq_len)
    return path
