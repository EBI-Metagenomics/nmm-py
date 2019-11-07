from typing import List, Tuple, Iterator, NamedTuple
from .._step import CStep
from .._path import CPath

Interval = NamedTuple("Interval", [("start", int), ("end", int)])


class Fragment:
    def __init__(
        self, seq: bytes, steps: List[CStep], homologous: bool, interval: Interval
    ):
        self._seq = seq
        self._homologous = homologous
        self._path = CPath(None)
        for step in steps:
            self._path.append(step.state.imm_state, step.seq_len)
        self._interval = interval

    @property
    def interval(self):
        return self._interval

    @property
    def sequence(self) -> bytes:
        return self._seq

    def items(self) -> Iterator[Tuple[bytes, CStep]]:
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
