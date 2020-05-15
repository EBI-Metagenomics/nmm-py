from typing import Generic, TypeVar

from .._cdata import CData
from .._ffi import ffi, lib
from ._hmm import HMM
from ._results import Results
from ._sequence import Sequence
from ._state import State

TState = TypeVar("TState", bound=State)


class DP(Generic[TState]):
    def __init__(self, imm_dp: CData, hmm: HMM[TState]):
        if imm_dp == ffi.NULL:
            raise RuntimeError("`imm_dp` is NULL.")
        self._imm_dp = imm_dp
        self._hmm = hmm

    def viterbi(self, seq: Sequence, window_length: int = 0) -> Results[TState]:
        from ._results import wrap_imm_results

        imm_seq = seq.imm_seq

        imm_results = lib.imm_dp_viterbi(self._imm_dp, imm_seq, window_length)
        if imm_results == ffi.NULL:
            raise RuntimeError("Could not run viterbi.")

        return wrap_imm_results(imm_results, seq, self._hmm.states())

    def __del__(self):
        if self._imm_dp != ffi.NULL:
            lib.imm_dp_destroy(self._imm_dp)
