from ._state import State, CState

from ._ffi import ffi, lib


class CStep:
    """
    Wrapper around the C implementation of path step.

    Parameters
    ----------
    imm_step : `<cdata 'struct imm_step *'>`.
    """

    def __init__(self, imm_step: ffi.CData):
        self.__cdata = imm_step

    @property
    def imm_step(self) -> ffi.CData:
        if self.__cdata == ffi.NULL:
            raise RuntimeError("`imm_step` is NULL.")
        return self.__cdata

    @property
    def state(self) -> CState:
        return CState(lib.imm_step_state(self.imm_step))

    @property
    def seq_len(self) -> int:
        return lib.imm_step_seq_len(self.imm_step)

    # def _set_imm_step(self, imm_step: ffi.CData):
    #     if self.__cdata != ffi.NULL:
    #         raise RuntimeError("`imm_step` is not NULL.")
    #     self.__cdata = imm_step

    def __str__(self) -> str:
        name = self.state.name.decode()
        return f"<{name},{self.seq_len}>"

    def __repr__(self) -> str:
        name = self.state.name.decode()
        return f"<{self.__class__.__name__}:{name},{self.seq_len}>"


class Step(CStep):
    """
    Path step.

    Parameters
    ----------
    state : State.
    seq_len : Sequence length.
    """

    def __init__(self, imm_step: ffi.CData, state: State, seq_len: int):
        super().__init__(imm_step)
        self._state = state
        self._seq_len = seq_len

    @property
    def state(self) -> State:
        return self._state

    @property
    def seq_len(self) -> int:
        return self._seq_len

    # def set_imm_step(self, imm_step: ffi.CData):
    #     self._set_imm_step(imm_step)
