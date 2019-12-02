from abc import ABC, abstractmethod
from typing import Iterator, NamedTuple, Sequence, Tuple

from .._path import PathBase
from .._step import StepBase

Interval = NamedTuple("Interval", [("start", int), ("end", int)])


class Fragment(ABC):
    def __init__(self, homologous: bool):
        self._homologous = homologous

    @property
    @abstractmethod
    def sequence(self) -> bytes:
        raise NotImplementedError()

    @abstractmethod
    def items(self) -> Iterator[Tuple[bytes, StepBase]]:
        raise NotImplementedError()

    @property
    def homologous(self) -> bool:
        return self._homologous


class SearchResult(ABC):
    def _create_fragments(self, path: PathBase):

        frag_start = frag_end = 0
        step_start = step_end = 0
        homologous = False

        for step_end, step in enumerate(path.steps()):

            change = not homologous and step.state.name.startswith(b"M")
            change = change or homologous and step.state.name.startswith(b"E")
            change = change or not homologous and step.state.name.startswith(b"T")

            if change:
                if frag_start < frag_end:
                    fragi = Interval(frag_start, frag_end)
                    stepi = Interval(step_start, step_end)
                    yield (fragi, stepi, homologous)

                frag_start = frag_end
                step_start = step_end
                homologous = not homologous

            frag_end += step.seq_len

    @property
    def sequence(self) -> bytes:
        return b"".join(frag.sequence for frag in self.fragments)

    @property
    @abstractmethod
    def fragments(self) -> Sequence[Fragment]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def intervals(self) -> Sequence[Interval]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def score(self) -> float:
        raise NotImplementedError()
