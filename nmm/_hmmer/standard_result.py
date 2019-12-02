from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval, SearchResult
from .standard_core import StandardPath, StandardStep


class StandardFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: StandardPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> bytes:
        return self._sequence

    def items(self) -> Iterator[Tuple[bytes, StandardStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class StandardSearchResult(SearchResult):
    def __init__(self, score: float, sequence: bytes, path: StandardPath):
        self._score = score

        self._fragments: List[StandardFragment] = []
        self._intervals: List[Interval] = []

        steps = list(path.steps())
        for fragi, stepi, homologous in self._create_fragments(path):
            spath = _create_path(steps[stepi.start : stepi.end])
            seq = sequence[fragi.start : fragi.end]
            frag = StandardFragment(seq, spath, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> Sequence[StandardFragment]:
        return self._fragments

    @property
    def interval(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def score(self) -> float:
        return self._score


def _create_path(steps: List[StandardStep]):
    path = StandardPath()
    for step in steps:
        path.append_standard_step(step.state, step.seq_len)
    return path
