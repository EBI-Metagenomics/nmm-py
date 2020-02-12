from typing import Sequence, Tuple, Union

from ..._ffi import lib
from ..._path import CPath
from ..._state import MuteState, NormalState
from .step import StandardStep


class StandardPath(CPath):
    """
    Path for standard profile.

    Parameters
    ----------
    steps : `Sequence[Tuple[Union[MuteState, NormalState], int]]`
        Steps.
    """

    def __init__(self, steps: Sequence[Tuple[Union[MuteState, NormalState], int]]):
        imm_path = lib.imm_path_create()
        self.__steps = [StandardStep(step[0], step[1]) for step in steps]
        for step in self.__steps:
            lib.imm_path_append(imm_path, step.imm_step)
        super().__init__(imm_path, self.__steps)

    def __getitem__(self, i) -> StandardStep:
        return self.__steps[i]

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
