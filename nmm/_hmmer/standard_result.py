from typing import Iterator, List, Sequence, Tuple

from .result import Fragment, Interval
from .standard_core import StandardPath, StandardStep


class StandardFragment(Fragment):
    def __init__(
        self,
        seq: bytes,
        steps: List[StandardStep],
        homologous: bool,
        interval: Interval,
    ):
        super().__init__(seq, homologous, interval)
        self._path = StandardPath()
        for step in steps:
            self._path.append_standard(step.state, step.seq_len)

    def items(self) -> Iterator[Tuple[bytes, StandardStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self._seq[start:end], step)
            start = end

    @property
    def homologous(self):
        return self._homologous

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class StandardSearchResult:
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
                    s = seq[frag_start:frag_end]
                    i = Interval(frag_start + 1, frag_end)
                    frag = StandardFragment(s, steps[idx_start:idx_end], False, i)
                    self._fragments.append(frag)
                homologous = True
                frag_start = frag_end
                idx_start = idx_end

            elif homologous and name.startswith(b"E"):
                if frag_start < frag_end:
                    s = seq[frag_start:frag_end]
                    i = Interval(frag_start + 1, frag_end)
                    frag = StandardFragment(s, steps[idx_start:idx_end], True, i)
                    self._fragments.append(frag)
                homologous = False
                frag_start = frag_end
                idx_start = idx_end

            frag_end += seq_len
            idx_end += 1

    @property
    def fragments(self) -> Sequence[StandardFragment]:
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
