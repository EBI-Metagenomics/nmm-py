from typing import List, Sequence, Dict

from .._cdata import CData
from ._state import State
from .._ffi import ffi, lib
from ._result import CResult
from ._sequence import Sequence as Seq


class CResults:
    """
    Wrapper around the C implementation of results.

    Parameters
    ----------
    imm_results : `<cdata 'struct imm_results *'>`.
        Results pointer.
    sequence : `Seq`
        Sequence.
    """

    def __init__(self, imm_results: CData, results: Sequence[CResult], sequence: Seq):
        if imm_results == ffi.NULL:
            raise RuntimeError("`imm_results` is NULL.")
        self._imm_results = imm_results
        self._results = list(results)
        self._sequence = sequence

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


def wrap_imm_results(imm_results: CData, sequence: Seq, states: Dict[CData, State]):
    from ._result import wrap_imm_result

    results: List[CResult] = []
    for i in range(lib.imm_results_size(imm_results)):
        imm_result = lib.imm_results_get(imm_results, i)
        results.append(wrap_imm_result(imm_result, sequence, states))

    return CResults(imm_results, results, sequence)
