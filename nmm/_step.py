from ._state import State

# from ._ffi import ffi, lib


class Step:
    def __init__(self, state: State, seq_len: int):
        self._state = state
        self._seq_len = seq_len

    @property
    def state(self) -> State:
        return self._state

    @property
    def seq_len(self) -> int:
        return self._seq_len
