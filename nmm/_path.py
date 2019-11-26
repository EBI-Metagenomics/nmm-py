from typing import List, Sequence, Tuple, Type, TypeVar, Optional

from ._ffi import ffi, lib
from ._step import CStep, Step
from ._state import State


class CPath:
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
            imm_path = lib.imm_path_create()
            if imm_path == ffi.NULL:
                raise RuntimeError("`imm_path_create` failed.")
        elif imm_path == ffi.NULL:
            raise RuntimeError("`cdata` cannot be NULL.")
        self.__cdata = imm_path

    @property
    def imm_path(self) -> ffi.CData:
        return self.__cdata

    def steps(self):
        step = lib.imm_path_first(self.__cdata)
        while step != ffi.NULL:
            yield CStep(step)
            step = lib.imm_path_next(self.__cdata, step)

    def append(self, state: ffi.CData, seq_len: int) -> ffi.CData:
        err: int = lib.imm_path_append(self.__cdata, state, seq_len)
        if err != 0:
            raise RuntimeError("Could not add step.")
        return lib.imm_path_last(self.__cdata)

    def __del__(self):
        if self.__cdata != ffi.NULL:
            lib.imm_path_destroy(self.__cdata)

    def __repr__(self):
        step_repr = []
        step = lib.imm_path_first(self.__cdata)
        while step != ffi.NULL:
            name = ffi.string(lib.imm_state_get_name(lib.imm_step_state(step)))
            seq_len = lib.imm_step_seq_len(step)
            step_repr += [f"<{name.decode()}:{seq_len}>"]
            step = lib.imm_path_next(self.__cdata, step)

        msg = "".join(step_repr)
        return f"<{self.__class__.__name__}:{msg}>"


T = TypeVar("T", bound="Path")


class Path(CPath):
    """
    Path of steps through a Markov model.

    Each step represents a state and an emitted sequence length.
    """

    def __init__(self):
        super().__init__()
        self._steps: List[Step] = []

    @classmethod
    def create(cls: Type[T], steps: Sequence[Tuple[State, int]]) -> T:
        path = cls()
        for state, seq_len in steps:
            path.append(state, seq_len)
        return path

    def append(self, state: State, seq_len: int) -> ffi.CData:
        imm_step = super().append(state.imm_state, seq_len)
        step = Step(imm_step, state, seq_len)
        self._steps.append(step)
        return step
