from typing import List, Sequence, Tuple

from ._ffi import ffi, lib
from ._state import CState
from ._step import CStep


class CPath:
    """
    Wrapper around the C implementation of path.

    Parameters
    ----------
    imm_path : `<cdata 'struct imm_path *'>`
        Path pointer.
    """

    def __init__(self, imm_path: ffi.CData):
        if imm_path == ffi.NULL:
            raise RuntimeError("`imm_path` is NULL.")
        self._imm_path = imm_path

        self._steps: List[CStep] = []
        step = lib.imm_path_first(imm_path)
        while step != ffi.NULL:
            self._steps.append(CStep(step))
            step = lib.imm_path_next(imm_path, step)

    @property
    def imm_path(self) -> ffi.CData:
        return self._imm_path

    def __len__(self) -> int:
        return len(self._steps)

    def __getitem__(self, i) -> CStep:
        return self._steps[i]

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
    def __init__(self, steps: Sequence[Tuple[CState, int]] = []):
        imm_path = lib.imm_path_create()

        for step in steps:
            if lib.imm_path_append(imm_path, step[0].imm_state, step[1]) != 0:
                lib.imm_path_destroy(imm_path)
                raise RuntimeError("Could not add step.")

        super().__init__(imm_path)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
