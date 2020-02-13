from ._state import CState
from .._ffi import ffi, lib


class CStep:
    """
    Wrapper around the C implementation of a path step.

    A step is composed of a state and an emitted sequence length. The user should not need to
    directly call the constructor of this class but instead use the methods from the `Path` class.

    Parameters
    ----------
    imm_step : `<cdata 'struct imm_step *'>`
        Step pointer.
    state : `CState`
        State.
    """

    def __init__(self, imm_step: ffi.CData, state: CState):
        if imm_step == ffi.NULL:
            raise RuntimeError("`imm_step` is NULL.")
        self._imm_step = imm_step
        self.__state = state

    @property
    def imm_step(self) -> ffi.CData:
        return self._imm_step

    @property
    def state(self) -> CState:
        return self.__state

    @property
    def seq_len(self) -> int:
        return lib.imm_step_seq_len(self.imm_step)

    def __del__(self):
        if self._imm_step != ffi.NULL:
            lib.imm_step_destroy(self._imm_step)

    def __str__(self) -> str:
        state = lib.imm_step_state(self._imm_step)
        name: str = ffi.string(lib.imm_state_get_name(state)).decode()
        return f"<{name},{self.seq_len}>"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Step(CStep):
    """
    Path step.

    Parameters
    ----------
    state : `CState`
        State.
    seq_len : `int`
        Sequence length.
    """

    def __init__(self, state: CState, seq_len: int):
        imm_step = lib.imm_step_create(state.imm_state, seq_len)
        if imm_step == ffi.NULL:
            raise RuntimeError("Could not create step.")
        super().__init__(imm_step, state)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"