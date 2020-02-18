from __future__ import annotations

from typing import (
    Dict,
    Generic,
    Iterable,
    Iterator,
    List,
    Sequence,
    Tuple,
    Type,
    TypeVar,
)

from .._ffi import ffi, lib
from ._state import CState
from ._step import CStep

S = TypeVar("S", bound=CState)
T = TypeVar("T", bound=CStep)


class CPath(Generic[T]):
    """
    Wrapper around the C implementation of path.

    Parameters
    ----------
    imm_path : `<cdata 'struct imm_path *'>`
        Path pointer.
    steps : `Iterable[T]`
        Steps.
    """

    def __init__(self, imm_path: ffi.CData, steps: Iterable[T]):
        if imm_path == ffi.NULL:
            raise RuntimeError("`imm_path` is NULL.")
        self._imm_path = imm_path
        self._steps = list(steps)

    @property
    def imm_path(self) -> ffi.CData:
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


class Path(CPath[T]):
    """
    Path.

    Parameters
    ----------
    steps : `Iterable[T]`
        Steps.
    """

    def __init__(self, steps: Iterable[T]):
        imm_path = lib.imm_path_create()
        for step in steps:
            lib.imm_path_append(imm_path, step.imm_step)
        super().__init__(imm_path, list(steps))

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


def wrap_imm_path(imm_path: ffi.CData, states: Dict[ffi.CData, CState]) -> CPath:
    steps: List[CStep] = []
    imm_step = lib.imm_path_first(imm_path)
    while imm_step != ffi.NULL:
        imm_state = lib.imm_step_state(imm_step)
        steps.append(CStep(imm_step, states[imm_state]))
        imm_step = lib.imm_path_next(imm_path, imm_step)

    return CPath(imm_path, steps)


def create_imm_path(steps: Sequence[CStep]) -> ffi.CData:
    imm_path = lib.imm_path_create()
    for step in steps:
        lib.imm_path_append(imm_path, step.imm_step)
    if imm_path == ffi.NULL:
        raise RuntimeError("Could not create path.")
    return imm_path
