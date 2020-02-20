from typing import Generic, Iterator, TypeVar

from .._interval import Interval
from ._alphabet import Alphabet
from ._path import Path
from ._sequence import SequenceABC
from ._state import State
from ._step import Step

TAlphabet = TypeVar("TAlphabet", bound=Alphabet)
TState = TypeVar("TState", bound=State)


class FragStep(Generic[TAlphabet, TState]):
    def __init__(self, sequence: SequenceABC[TAlphabet], step: Step[TState]):
        self._sequence = sequence
        self._step = step

    @property
    def sequence(self) -> SequenceABC[TAlphabet]:
        return self._sequence

    @property
    def step(self) -> Step[TState]:
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
        self, sequence: SequenceABC[TAlphabet], path: Path[Step[TState]],
    ):
        self._sequence = sequence
        self._path = path

    @property
    def sequence(self) -> SequenceABC[TAlphabet]:
        return self._sequence

    @property
    def path(self) -> Path[Step[TState]]:
        return self._path

    def __iter__(self) -> Iterator[FragStep]:
        start = end = 0
        for step in self._path:
            end += step.seq_len
            yield FragStep(self._sequence[Interval(start, end)], step)
            start = end

    def __str__(self) -> str:
        return ",".join(str(seq_step) for seq_step in self)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
