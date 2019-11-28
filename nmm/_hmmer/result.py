from typing import Iterator, NamedTuple, Sequence, Tuple

from .._step import Step

Interval = NamedTuple("Interval", [("start", int), ("end", int)])


class Fragment:
    def __init__(self, seq: bytes, homologous: bool, interval: Interval):
        self._seq = seq
        self._homologous = homologous
        self._interval = interval

    @property
    def interval(self) -> Interval:
        return self._interval

    @property
    def sequence(self) -> bytes:
        return self._seq

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
        return self._fragments

    @property
    def score(self) -> float:
        raise NotImplementedError()
