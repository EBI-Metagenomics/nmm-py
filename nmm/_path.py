from typing import List, Sequence, Union

from ._ffi import ffi, lib
from ._step import CStep, Step


class CPath:
    def __init__(self, cdata: Union[ffi.CData, None]):
        if cdata is None:
            cdata = lib.imm_path_create()
            if cdata == ffi.NULL:
                raise RuntimeError("`imm_path_create` failed.")
        elif cdata == ffi.NULL:
            raise RuntimeError("`cdata` cannot be NULL.")
        self.__cdata = cdata

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
        if self._cdata != ffi.NULL:
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


class Path(CPath):
    def __init__(self, steps: Sequence[Step]):
        super().__init__(None)

        self._steps: List[CStep] = []
        for step in steps:
            imm_step = self.append(step.state.imm_state, step.seq_len)
            step.set_imm_step(imm_step)
            self._steps.append(step)
