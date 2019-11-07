from typing import List, Union

from .._path import CPath
from .fragment import Fragment, Interval


class Result:
    def __init__(self, score: float, seq: bytes, path: CPath):
        self._score = score

        self._fragments: List[Union[Fragment]] = []

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
                    frag = Fragment(s, steps[idx_start:idx_end], False, i)
                    self._fragments.append(frag)
                homologous = True
                frag_start = frag_end
                idx_start = idx_end

            elif homologous and name.startswith(b"E"):
                if frag_start < frag_end:
                    s = seq[frag_start:frag_end]
                    i = Interval(frag_start + 1, frag_end)
                    frag = Fragment(s, steps[idx_start:idx_end], True, i)
                    self._fragments.append(frag)
                homologous = False
                frag_start = frag_end
                idx_start = idx_end

            frag_end += seq_len
            idx_end += 1

    @property
    def fragments(self):
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
