from __future__ import annotations

from typing import Dict, Generic, Iterable, Iterator, List, Type, TypeVar

from .._cdata import CData
from .._ffi import ffi, lib
from ._state import State
from ._step import Step

T = TypeVar("T", bound=Step)


class Path(Generic[T]):
    """
    Path.

    Parameters
    ----------
    imm_path
        Path pointer.
    steps
        Steps.
    """

    def __init__(self, imm_path: CData, steps: Iterable[T]):
        if imm_path == ffi.NULL:
            raise RuntimeError("`imm_path` is NULL.")
        self._imm_path = imm_path
        self._steps = list(steps)

    @classmethod
    def create(cls: Type[Path], steps: Iterable[T]) -> Path:
        """
        Create path.

        Parameters
        ----------
        steps
            Steps.
        """
        imm_path = lib.imm_path_create()
        for step in steps:
            lib.imm_path_append(imm_path, step.imm_step)
        return cls(imm_path, list(steps))

    @property
    def imm_path(self) -> CData:
        return self._imm_path

    def __len__(self) -> int:
        return len(self._steps)

    def __getitem__(self, i: int) -> T:
        return self._steps[i]

    def __iter__(self) -> Iterator[T]:
        for i in range(len(self)):
            yield self[i]

    def __del__(self):
        if self._imm_path != ffi.NULL:
            lib.imm_path_free(self._imm_path)

    def __str__(self) -> str:
        return ",".join([str(s) for s in self])

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


def wrap_imm_path(imm_path: CData, states: Dict[CData, State]) -> Path:
    steps: List[Step] = []
    imm_step = lib.imm_path_first(imm_path)
    while imm_step != ffi.NULL:
        imm_state = lib.imm_step_state(imm_step)
        steps.append(Step(imm_step, states[imm_state]))
        imm_step = lib.imm_path_next(imm_path, imm_step)

    return Path(imm_path, steps)
