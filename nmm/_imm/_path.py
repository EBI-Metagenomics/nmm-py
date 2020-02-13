from typing import Dict, List, Sequence, Tuple

from .._ffi import ffi, lib
from ._state import CState
from ._step import CStep, Step


class CPath:
    """
    Wrapper around the C implementation of path.

    Parameters
    ----------
    imm_path : `<cdata 'struct imm_path *'>`
        Path pointer.
    steps : `Sequence[CStep]`
        List of steps.
    """

    def __init__(self, imm_path: ffi.CData, steps: Sequence[CStep]):
        if imm_path == ffi.NULL:
            raise RuntimeError("`imm_path` is NULL.")
        self._imm_path = imm_path
        self.__steps = list(steps)

    @property
    def imm_path(self) -> ffi.CData:
        return self._imm_path

    def __len__(self) -> int:
        return len(self.__steps)

    def __getitem__(self, i) -> CStep:
        return self.__steps[i]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __del__(self):
        if self._imm_path != ffi.NULL:
            lib.imm_path_free(self._imm_path)

    def __str__(self) -> str:
        return ",".join([str(s) for s in self])

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class Path(CPath):
    """
    Path.

    Parameters
    ----------
    steps : `Sequence[Tuple[CState, int]]`
        Steps.
    """

    def __init__(self, steps: Sequence[Tuple[CState, int]]):
        imm_path = lib.imm_path_create()
        self.__steps = [Step(step[0], step[1]) for step in steps]
        for step in self.__steps:
            lib.imm_path_append(imm_path, step.imm_step)
        super().__init__(imm_path, self.__steps)

    def __getitem__(self, i) -> CStep:
        return self.__steps[i]

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
