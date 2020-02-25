from typing import Mapping, TypeVar, Generic

from .._cdata import CData
from .._ffi import ffi, lib
from ._path import Path, wrap_imm_path
from ._step import Step
from ._sequence import Sequence, SequenceABC, SubSequence
from ._state import State


TState = TypeVar("TState", bound=State)


class Result(Generic[TState]):
    """
    Result.

    Parameters
    ----------
    imm_result
        Result pointer.
    path
        Path.
    sequence
        Sequence.
    """

    def __init__(
        self, imm_result: CData, path: Path[Step[TState]], sequence: SequenceABC
    ):
        if imm_result == ffi.NULL:
            raise RuntimeError("`imm_result` is NULL.")
        self._imm_result = imm_result
        self._path = path
        self._sequence = sequence

    @property
    def loglikelihood(self) -> float:
        return lib.imm_result_loglik(self._imm_result)

    @property
    def path(self) -> Path[Step[TState]]:
        return self._path

    @property
    def sequence(self) -> SequenceABC:
        return self._sequence

    def __del__(self):
        if self._imm_result != ffi.NULL:
            lib.imm_result_free(self._imm_result)

    def __repr__(self) -> str:
        return str(self.loglikelihood)


def wrap_imm_result(
    imm_result: CData, sequence: Sequence, states: Mapping[CData, TState]
):
    path = wrap_imm_path(lib.imm_result_path(imm_result), states)
    imm_subseq = lib.imm_result_subseq(imm_result)
    return Result(imm_result, path, SubSequence(imm_subseq, sequence))
