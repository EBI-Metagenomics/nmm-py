from typing import Union

from ..._ffi import ffi, lib
from ..._state import MuteState, NormalState
from ..._step import CStep


class StandardStep(CStep):
    """
    Path step for the standard profile.

    Parameters
    ----------
    state : `CState`
        State.
    seq_len : `int`
        Sequence length.
    """

    def __init__(self, state: Union[MuteState, NormalState], seq_len: int):
        imm_step = lib.imm_step_create(state.imm_state, seq_len)
        if imm_step == ffi.NULL:
            raise RuntimeError("Could not create step.")
        super().__init__(imm_step, state)
        self._state = state

    @property
    def state(self) -> Union[MuteState, NormalState]:
        return self._state

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
