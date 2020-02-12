from ._ffi import ffi, lib
from ._sequence import CSequence
from ._interval import Interval


class CSubSequence:
    """
    Wrapper around the C implementation of subsequence.

    Parameters
    ----------
    imm_subseq : `<cdata 'struct imm_subseq *'>`.
        Subsequence pointer.
    sequence : `CSequence`
        Sequence.
    """

    def __init__(self, imm_subseq: ffi.CData, sequence: CSequence):
        if ffi.getctype(ffi.typeof(imm_subseq)) != "struct imm_subseq":
            raise TypeError("Wrong `imm_subseq` type.")
        self._imm_subseq = imm_subseq
        self._sequence = sequence

    @property
    def imm_subseq(self) -> ffi.CData:
        return self._imm_subseq

    @property
    def start(self) -> int:
        return lib.imm_subseq_start(ffi.addressof(self._imm_subseq))

    @property
    def length(self) -> int:
        return lib.imm_subseq_length(ffi.addressof(self._imm_subseq))

    @property
    def symbols(self) -> bytes:
        imm_seq = lib.imm_subseq_cast(ffi.addressof(self._imm_subseq))
        return ffi.string(lib.imm_seq_string(imm_seq), lib.imm_seq_length(imm_seq))

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SubSequence(CSubSequence):
    """
    Subsequence of symbols of a given sequence.

    Parameters
    ----------
    sequence : `CSequence`
        Sequence.
    interval : `Interval`
        Interval.
    """

    def __init_(self, sequence: CSequence, interval: Interval):
        length = interval.stop - interval.start
        if interval.start < 0 or interval.start > interval.stop or length > self.length:
            raise ValueError("Out-of-range interval.")

        imm_subseq = lib.imm_subseq_slice(sequence.imm_seq, interval.start, length)
        super().__init__(imm_subseq, sequence)
