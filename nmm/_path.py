from typing import Sequence, Union, List
from ._step import Step

from ._ffi import ffi, lib


class CPath:
    def __init__(self, cdata: Union[ffi.CData, None]):
        # self._own = True
        self._cdata = ffi.NULL
        if cdata is None:
            self._cdata = lib.imm_path_create()
            if self._cdata == ffi.NULL:
                raise RuntimeError("`imm_path_create` failed.")
        else:
            self._cdata = cdata

    @property
    def cdata(self) -> ffi.CData:
        # if not self._own:
        # raise RuntimeError("I do not own this data anymore.")
        return self._cdata

    # def borrow_cdata(self) -> ffi.CData:
    #     cdata = self.cdata
    #     self._own = False
    #     return cdata

    def _append(self, state: ffi.CData, seq_len: int):
        err: int = lib.imm_path_append(self._cdata, state, seq_len)
        if err != 0:
            raise RuntimeError("Could not add step.")

    def __del__(self):
        # if self._own and self._cdata != ffi.NULL:
        if self._cdata != ffi.NULL:
            lib.imm_path_destroy(self._cdata)

    def __repr__(self):
        step = lib.imm_path_first(self._path)
        steps = []
        while step != ffi.NULL:
            name = lib.imm_state_get_name(lib.imm_step_state(step))
            state_name = ffi.string(name).decode()
            seq_len = lib.imm_step_seq_len(step)
            steps += [f"<{state_name}:{seq_len}>"]
            step = lib.imm_path_next(self._path, step)

        step_msg = "".join(steps)
        return f"<{self.__class__.__name__}:{step_msg}>"


class Path(CPath):
    def __init__(self, steps: Sequence[Step]):
        super().__init__(None)
        self._steps: List[Step] = []
        for step in steps:
            if step.seq_len < 0:
                raise ValueError("Sequence length cannot be negative.")
            self._append(step.state.imm_state, step.seq_len)
            self._steps.append(step)

    # def steps(self):
    #     step = lib.imm_path_first(self._path)
    #     while step != ffi.NULL:
    #         cname = lib.imm_state_get_name(lib.imm_step_state(step))
    #         name = ffi.string(cname).decode()
    #         seq_len = lib.imm_step_seq_len(step)

    #         step = lib.imm_path_next(self._path, step)
