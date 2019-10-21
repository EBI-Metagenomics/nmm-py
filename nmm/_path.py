from ._state import State

from ._ffi import ffi, lib


class Path:
    def __init__(self, path: list):
        self._path = lib.imm_path_create()
        for step in path:
            self._add(step[0].cdata, step[1])

    @property
    def cdata(self):
        return self._path

    def _add(self, state: State, seq_len: int):
        if seq_len < 0:
            raise ValueError("Sequence length cannot be negative.")
        lib.imm_path_add(self._path, state, seq_len)

    def __del__(self):
        if self._path != ffi.NULL:
            lib.imm_path_destroy(self._path)
