from abc import ABC, abstractmethod
from typing import Iterator, Optional, Sequence, Tuple, Type, TypeVar

from ._ffi import ffi, lib
from ._state import CState
from ._step import CStep, StepBase


class PathBase(ABC):
    @abstractmethod
    def steps(self) -> Iterator[StepBase]:
        raise NotImplementedError()


T = TypeVar("T", bound="CPath")


class CPath(PathBase):
    """
    Wrapper around the C implementation of path.

    Parameters
    ----------
    imm_path : Optional[ffi.CData]
        `<cdata 'struct imm_path *'>` or `None`. Passing `None` will create a new path at the
        underlying library level.
    """

    def __init__(self, imm_path: Optional[ffi.CData] = None):
        if imm_path is None:
            imm_path = _create_imm_path()
        self._imm_path = imm_path

    @classmethod
    def create_cpath(cls: Type[T], steps: Sequence[Tuple[CState, int]] = []) -> T:
        cpath = cls()
        for step in steps:
            cpath.append_cstep(step[0], step[1])
        return cpath

    @property
    def imm_path(self) -> ffi.CData:
        return self._imm_path

    def steps(self) -> Iterator[CStep]:
        step = lib.imm_path_first(self._imm_path)
        while step != ffi.NULL:
            yield CStep(step)
            step = lib.imm_path_next(self._imm_path, step)

    def append_cstep(self, state: CState, seq_len: int) -> CStep:

        err: int = lib.imm_path_append(self._imm_path, state.imm_state, seq_len)
        if err != 0:
            raise RuntimeError("Could not add step.")

        imm_step = lib.imm_path_last(self._imm_path)
        return CStep(imm_step)

    def __del__(self):
        if self._imm_path != ffi.NULL:
            lib.imm_path_destroy(self._imm_path)

    def __repr__(self):
        step_repr = [f"<{s.state.name.decode()}:{s.seq_len}>" for s in self.steps()]
        msg = "".join(step_repr)
        return f"<{self.__class__.__name__}:{msg}>"


def _create_imm_path():
    imm_path = lib.imm_path_create()
    if imm_path == ffi.NULL:
        raise RuntimeError("`imm_path_create` failed.")
    return imm_path
