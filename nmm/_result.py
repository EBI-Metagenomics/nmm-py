from ._path import CPath
from ._sequence import CSequence
from ._ffi import ffi, lib


class CResult:
    def __init__(self, imm_result: ffi.CData):
        if imm_result == ffi.NULL:
            raise RuntimeError("`imm_result` is NULL.")
        self._imm_result = imm_result

        self._path = CPath(lib.imm_path_clone(lib.imm_result_path(imm_result)))
        imm_seq = lib.imm_seq_clone(lib.imm_result_sequence(imm_result))
        self._sequence = CSequence(imm_seq)

    @property
    def loglikelihood(self) -> float:
        return lib.imm_result_loglik(self._imm_result)

    @property
    def path(self) -> CPath:
        return self._path

    @property
    def sequence(self) -> CSequence:
        return self._sequence

    def __del__(self):
        if self._imm_result != ffi.NULL:
            lib.imm_result_free(self._imm_result)

    def __repr__(self) -> str:
        return str(self.loglikelihood)
