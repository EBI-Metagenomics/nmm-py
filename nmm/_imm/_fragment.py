from typing import Generic, Iterator, TypeVar

from .._interval import Interval
from ._alphabet import CAlphabet
from ._path import CPath
from ._sequence import SequenceABC
from ._state import CState
from ._step import CStep

TAlphabet = TypeVar("TAlphabet", bound=CAlphabet)
TState = TypeVar("TState", bound=CState)


class FragStep(Generic[TAlphabet, TState]):
    def __init__(self, sequence: SequenceABC[TAlphabet], step: CStep[TState]):
        self._sequence = sequence
        self._step = step

    @property
    def sequence(self) -> SequenceABC[TAlphabet]:
        return self._sequence

    @property
    def step(self) -> CStep[TState]:
        return self._step

    def __str__(self) -> str:
        return f"{str(self.sequence), str(self.step)}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Fragment(Generic[TAlphabet, TState]):
    """
    Fragment of the standard profile.

    Parameters
    ----------
    sequence
        Sequence.
    path
        Path of the standard profile.
    """

    def __init__(
        self, sequence: SequenceABC[TAlphabet], path: CPath[CStep[TState]],
    ):
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> SequenceABC[TAlphabet]:
        return self._sequence

    def __iter__(self) -> Iterator[FragStep]:
        start = end = 0
        for step in self._path:
            end += step.seq_len
            yield FragStep(self._sequence.slice(Interval(start, end)), step)
            start = end

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
