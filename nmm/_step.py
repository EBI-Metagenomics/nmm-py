from ._state import State, CState

from ._ffi import ffi, lib


class CStep:
    def __init__(self, cdata: ffi.CData):
        self.__cdata = cdata

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

    def _set_imm_step(self, imm_step: ffi.CData):
        if self.__cdata != ffi.NULL:
            raise RuntimeError("`imm_step` is not NULL.")
        self.__cdata = imm_step


class Step(CStep):
    def __init__(self, state: State, seq_len: int):
        self._state = state
        self._seq_len = seq_len
        super().__init__(ffi.NULL)

    @property
    def state(self) -> State:
        return self._state

    @property
    def seq_len(self) -> int:
        return self._seq_len

    def set_imm_step(self, imm_step: ffi.CData):
        self._set_imm_step(imm_step)
