from typing import Iterator, NamedTuple, Sequence, Tuple

from .._step import Step

Interval = NamedTuple("Interval", [("start", int), ("end", int)])


class Fragment:
    def __init__(self, seq: bytes, interval: Interval, homologous: bool):
        self._seq = seq
        self._interval = interval
        self._homologous = homologous

    @property
    def interval(self) -> Interval:
        return self._interval

    @property
    def sequence(self) -> bytes:
        return self._seq[self._interval.start : self._interval.end]

    def items(self) -> Iterator[Tuple[bytes, Step]]:
        raise NotImplementedError()

    @property
    def homologous(self) -> bool:
        return self._homologous

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class SearchResult:
    @property
    def fragments(self) -> Sequence[Fragment]:
        raise NotImplementedError()

    @property
    def score(self) -> float:
        raise NotImplementedError()
