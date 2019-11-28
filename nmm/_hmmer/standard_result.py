from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval, SearchResult
from .standard_core import StandardPath, StandardStep


class StandardFragment(Fragment):
    def __init__(
        self, seq: bytes, interval: Interval, path: StandardPath, homologous: bool,
    ):
        super().__init__(seq, interval, homologous)
        self._path = path

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
    def __init__(self, score: float, seq: bytes, path: StandardPath):
        self._score = score

        self._fragments: List[StandardFragment] = []

        steps = list(path.steps())
        for fragi, stepi, homologous in self._create_fragments(path):
            spath = _create_path(steps[stepi.start : stepi.end])
            self._fragments.append(StandardFragment(seq, fragi, spath, homologous))

    @property
    def fragments(self) -> Sequence[StandardFragment]:
        return self._fragments

    @property
    def score(self) -> float:
        return self._score


def _create_path(steps: List[StandardStep]):
    path = StandardPath()
    for step in steps:
        path.append_standard_step(step.state, step.seq_len)
    return path
