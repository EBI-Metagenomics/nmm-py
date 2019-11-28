from typing import Iterator, NamedTuple, Sequence, Tuple

from .._step import Step
from .._path import Path

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
    def _create_fragments(self, path: Path):

        frag_start = frag_end = 0
        step_start = step_end = 0
        homologous = False

        for step_end, step in enumerate(path.steps()):

            change = not homologous and step.state.name.startswith(b"M")
            change = change or homologous and step.state.name.startswith(b"E")

            if change:
                if frag_start < frag_end:
                    fragi = Interval(frag_start, frag_end)
                    stepi = Interval(step_start, step_end)
                    yield (fragi, stepi, homologous)

                homologous = not homologous
                frag_start = frag_end
                step_start = step_end

            frag_end += step.seq_len

    @property
    def fragments(self) -> Sequence[Fragment]:
        raise NotImplementedError()

    @property
    def score(self) -> float:
        raise NotImplementedError()
