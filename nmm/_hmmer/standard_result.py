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


def _create_path(steps: List[StandardStep]):
    path = StandardPath()
    for step in steps:
        path.append_standard_step(step.state, step.seq_len)
    return path


class StandardSearchResult(SearchResult):
    def __init__(self, score: float, seq: bytes, path: StandardPath):
        self._score = score

        self._fragments: List[StandardFragment] = []

        frag_start = frag_end = 0
        idx_start = idx_end = 0
        homologous = False

        steps = list(path.steps())
        for step in steps:
            name = step.state.name
            seq_len = step.seq_len

            if not homologous and name.startswith(b"M"):
                if frag_start < frag_end:
                    # s = seq[frag_start:frag_end]
                    i = Interval(frag_start, frag_end)
                    spath = _create_path(steps[idx_start:idx_end])
                    frag = StandardFragment(seq, i, spath, False)
                    self._fragments.append(frag)
                homologous = True
                frag_start = frag_end
                idx_start = idx_end

            elif homologous and name.startswith(b"E"):
                if frag_start < frag_end:
                    # s = seq[frag_start:frag_end]
                    i = Interval(frag_start, frag_end)
                    spath = _create_path(steps[idx_start:idx_end])
                    frag = StandardFragment(seq, i, spath, True)
                    self._fragments.append(frag)
                homologous = False
                frag_start = frag_end
                idx_start = idx_end

            frag_end += seq_len
            idx_end += 1

    # def _append_fragment(self, seq: bytes, homologous: bool, ):
    #     pass

    @property
    def fragments(self) -> Sequence[StandardFragment]:
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
