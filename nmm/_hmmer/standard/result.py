from typing import List, Sequence

from ..result import SearchResult
from .fragment import StandardFragment
from .path import StandardPath
from ..._sequence import CSequence
from ..._interval import Interval


class StandardSearchResult(SearchResult):
    def __init__(self, loglik: float, sequence: CSequence, path: StandardPath):
        self._loglik = loglik
        self._fragments: List[StandardFragment] = []
        self._intervals: List[Interval] = []

        steps = list(path)
        for fragi, stepi, homologous in self._create_fragments(path):
            # fragment_path = _create_path(steps[stepi.start : stepi.stop])
            substeps = steps[stepi.start : stepi.stop]
            fragment_path = StandardPath([(s.state, s.seq_len) for s in substeps])
            seq = sequence[fragi.start : fragi.stop]
            frag = StandardFragment(seq, fragment_path, homologous)
            self._fragments.append(frag)
            self._intervals.append(fragi)

    @property
    def fragments(self) -> Sequence[StandardFragment]:
        return self._fragments

    @property
    def intervals(self) -> Sequence[Interval]:
        return self._intervals

    @property
    def loglikelihood(self) -> float:
        return self._loglik
