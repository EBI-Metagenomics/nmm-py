from typing import Dict

from .._ffi import ffi, lib
from ._path import CPath, wrap_imm_path
from ._sequence import CSequence, CSubSequence
from ._state import CState


class CResult:
    def __init__(self, imm_result: ffi.CData, path: CPath, subseq: CSubSequence):
        if imm_result == ffi.NULL:
            raise RuntimeError("`imm_result` is NULL.")
        self._imm_result = imm_result
        self._path = path
        self._subseq = subseq

    @property
    def loglikelihood(self) -> float:
        return lib.imm_result_loglik(self._imm_result)

    @property
    def path(self) -> CPath:
        return self._path

    @property
    def subseq(self) -> CSubSequence:
        return self._subseq

    def __del__(self):
        if self._imm_result != ffi.NULL:
            lib.imm_result_free(self._imm_result)

    def __repr__(self) -> str:
        return str(self.loglikelihood)


def wrap_imm_result(
    imm_result: ffi.CData, sequence: CSequence, states: Dict[ffi.CData, CState]
):
    path = wrap_imm_path(lib.imm_result_path(imm_result), states)
    imm_subseq = lib.imm_result_subseq(imm_result)
    subseq = CSubSequence(imm_subseq, sequence)
    return CResult(imm_result, path, subseq)
