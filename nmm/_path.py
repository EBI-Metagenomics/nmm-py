from typing import Sequence, Tuple
from ._state import State

from ._ffi import ffi, lib


class Path:
    def __init__(self, path: Sequence[Tuple[State, int]]):
        self._path = lib.imm_path_create()
        for step in path:
            self._append(step[0].cdata, step[1])

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
