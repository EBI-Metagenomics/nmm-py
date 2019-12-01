from abc import ABC, abstractmethod
from typing import List, Sequence, Tuple, Type, TypeVar, Optional, Iterator

from ._ffi import ffi, lib
from ._step import CStep, StepBase
from ._state import StateBase, CState


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

    # def __init__(self, imm_path: Optional[ffi.CData] = None):
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

    # def _append_imm_step(self, state: ffi.CData, seq_len: int) -> ffi.CData:
    #     err: int = lib.imm_path_append(self._imm_path, state, seq_len)
    #     if err != 0:
    #         raise RuntimeError("Could not add step.")
    #     return lib.imm_path_last(self._imm_path)

    def __del__(self):
        if self._imm_path != ffi.NULL:
            lib.imm_path_destroy(self._imm_path)

    def __repr__(self):
        step_repr = []
        step = lib.imm_path_first(self._imm_path)
        while step != ffi.NULL:
            name = ffi.string(lib.imm_state_get_name(lib.imm_step_state(step)))
            seq_len = lib.imm_step_seq_len(step)
            step_repr += [f"<{name.decode()}:{seq_len}>"]
            step = lib.imm_path_next(self._imm_path, step)

        msg = "".join(step_repr)
        return f"<{self.__class__.__name__}:{msg}>"


# class CPath(PathBase):
#     """
#     Wrapper around the C implementation of path.

#     Parameters
#     ----------
#     imm_path : Optional[ffi.CData]
#         `<cdata 'struct imm_path *'>` or `None`. Passing `None` will create a new path at the
#         underlying library level.
#     """

#     def __init__(self, imm_path: Optional[ffi.CData] = None):
#         if imm_path is None:
#             imm_path = lib.imm_path_create()
#             if imm_path == ffi.NULL:
#                 raise RuntimeError("`imm_path_create` failed.")
#         elif imm_path == ffi.NULL:
#             raise RuntimeError("`cdata` cannot be NULL.")
#         self.__imm_path = imm_path

#     @property
#     def imm_path(self) -> ffi.CData:
#         return self.__imm_path

#     def steps(self) -> Iterator[CStep]:
#         step = lib.imm_path_first(self.__imm_path)
#         while step != ffi.NULL:
#             yield CStep(step)
#             step = lib.imm_path_next(self.__imm_path, step)

#     def _append_imm_step(self, state: ffi.CData, seq_len: int) -> ffi.CData:
#         err: int = lib.imm_path_append(self.__imm_path, state, seq_len)
#         if err != 0:
#             raise RuntimeError("Could not add step.")
#         return lib.imm_path_last(self.__imm_path)

#     def __del__(self):
#         if self.__imm_path != ffi.NULL:
#             lib.imm_path_destroy(self.__imm_path)

#     def __repr__(self):
#         step_repr = []
#         step = lib.imm_path_first(self.__imm_path)
#         while step != ffi.NULL:
#             name = ffi.string(lib.imm_state_get_name(lib.imm_step_state(step)))
#             seq_len = lib.imm_step_seq_len(step)
#             step_repr += [f"<{name.decode()}:{seq_len}>"]
#             step = lib.imm_path_next(self.__imm_path, step)

#         msg = "".join(step_repr)
#         return f"<{self.__class__.__name__}:{msg}>"


# class Path(CPath):
#     """
#     Path of steps through a Markov model.

#     Each step represents a state and an emitted sequence length.
#     """

#     def __init__(self):
#         super().__init__()
#         self.__steps: List[Step] = []

#     @classmethod
#     def create(cls: Type[T], steps: Sequence[Tuple[CState, int]]) -> T:
#         path = cls()
#         for state, seq_len in steps:
#             path.append(state, seq_len)
#         return path

#     def append(self, state: CState, seq_len: int) -> ffi.CData:
#         imm_step = self._append_imm_step(state.imm_state, seq_len)
#         step = Step(imm_step, state)
#         self.__steps.append(step)
#         return step

#     def steps(self) -> Iterator[Step]:
#         return iter(self.__steps)


def _create_imm_path():
    imm_path = lib.imm_path_create()
    if imm_path == ffi.NULL:
        raise RuntimeError("`imm_path_create` failed.")
    return imm_path
