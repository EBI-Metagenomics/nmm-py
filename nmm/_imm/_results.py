from typing import Mapping, Iterable, List, TypeVar, Generic

from .._cdata import CData
from .._ffi import ffi, lib
from ._result import Result
from ._sequence import Sequence, SequenceABC
from ._state import State

TState = TypeVar("TState", bound=State)


class Results(Generic[TState]):
    """
    Results.

    Parameters
    ----------
    imm_results
        Results pointer.
    results
        List of results.
    sequence
        Sequence.
    """

    def __init__(
        self,
        imm_results: CData,
        results: Iterable[Result[TState]],
        sequence: SequenceABC,
    ):
        if imm_results == ffi.NULL:
            raise RuntimeError("`imm_results` is NULL.")
        self._imm_results = imm_results
        self._results = list(results)
        self._sequence = sequence

    def __len__(self) -> int:
        return len(self._results)

    def __getitem__(self, i) -> Result[TState]:
        return self._results[i]

    def __iter__(self):
        for r in self._results:
            yield r

    def __del__(self):
        if self._imm_results != ffi.NULL:
            lib.imm_results_free(self._imm_results)

    def __repr__(self) -> str:
        return "[" + ",".join([str(r) for r in self]) + "]"


def wrap_imm_results(
    imm_results: CData, sequence: Sequence, states: Mapping[CData, TState]
):
    from ._result import wrap_imm_result

    results: List[Result[TState]] = []
    for i in range(lib.imm_results_size(imm_results)):
        imm_result = lib.imm_results_get(imm_results, i)
        results.append(wrap_imm_result(imm_result, sequence, states))

    return Results(imm_results, results, sequence)
