from typing import Sequence, Tuple, Union
from ._state import State

from ._ffi import ffi, lib


class CPath:
    def __init__(self):
        self._own = True
        self._cdata = ffi.NULL
        self._cdata = lib.imm_path_create()
        if self._cdata == ffi.NULL:
            raise RuntimeError("Could not create CPath.")

    @property
    def cdata(self) -> ffi.CData:
        if not self._own:
            raise RuntimeError("I do not own this data anymore.")
        return self._cdata

    def borrow_cdata(self) -> ffi.CData:
        cdata = self.cdata
        self._own = False
        return cdata

    def __del__(self):
        if self._own and self._cdata != ffi.NULL:
            lib.imm_path_destroy(self._cdata)


class Path:
    def __init__(self, steps: Union[Sequence[Tuple[State, int]], CPath]):
        self._path = ffi.NULL
        if isinstance(steps, CPath):
            self._path = steps.borrow_cdata()
        else:
            self._path = lib.imm_path_create()
            for step in steps:
                self._append(step[0].cdata, step[1])

    # def steps(self):
    #     step = lib.imm_path_first(self._path)
    #     while step != ffi.NULL:
    #         cname = lib.imm_state_get_name(lib.imm_step_state(step))
    #         name = ffi.string(cname).decode()
    #         seq_len = lib.imm_step_seq_len(step)

    #         step = lib.imm_path_next(self._path, step)

    @property
    def cdata(self):
        return self._path

    def _append(self, state: State, seq_len: int):
        if seq_len < 0:
            raise ValueError("Sequence length cannot be negative.")
        err: int = lib.imm_path_append(self._path, state, seq_len)
        if err != 0:
            raise RuntimeError("Could not add step.")

    def __del__(self):
        if self._path != ffi.NULL:
            lib.imm_path_destroy(self._path)

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

