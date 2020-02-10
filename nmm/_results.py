from typing import List

from ._ffi import ffi, lib
from ._result import CResult
from ._sequence import CSequence


class CResults:
    """
    Wrapper around the C implementation of results.

    Parameters
    ----------
    imm_results : `<cdata 'struct imm_results *'>`.
        Results pointer.
    sequence : `CSequence`
        Sequence.
    """

    def __init__(self, imm_results: ffi.CData, sequence: CSequence):
        if imm_results == ffi.NULL:
            raise RuntimeError("`imm_results` is NULL.")
        self._imm_results = imm_results
        self._sequence = sequence

        self._results: List[CResult] = []
        for i in range(lib.imm_results_size(self._imm_results)):
            self._results.append(CResult(lib.imm_results_get(self._imm_results, i)))

    def __len__(self) -> int:
        return len(self._results)

    def __getitem__(self, i) -> CResult:
        return self._results[i]

    def __iter__(self):
        for r in self._results:
            yield r

    def __del__(self):
        if self._imm_results != ffi.NULL:
            lib.imm_results_free(self._imm_results)

    def __repr__(self) -> str:
        return "[" + ",".join([str(r) for r in self]) + "]"
