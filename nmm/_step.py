from typing import Sequence, Tuple
from ._state import State

from ._ffi import ffi, lib


class Step:
    def __init__(self, state: State, seq_len: int):
        pass
