from __future__ import annotations

from typing import TypeVar, Generic
from abc import ABC, abstractmethod

from ._alphabet import CAlphabet
from .._ffi import ffi, lib
from .._interval import Interval


T = TypeVar("T", bound=CAlphabet)


class SequenceABC(ABC, Generic[T]):
    """
    Sequence of symbols.
    """

    @property
    @abstractmethod
    def length(self) -> int:
        raise NotImplementedError()

    @property
    @abstractmethod
    def symbols(self) -> bytes:
        raise NotImplementedError()

    @abstractmethod
    def slice(self, interval: Interval) -> SequenceABC[T]:
        del interval
        raise NotImplementedError()

    @property
    @abstractmethod
    def alphabet(self) -> T:
        raise NotImplementedError()


class CSequence(SequenceABC[T]):
    """
    Wrapper around the C implementation of sequence.

    Parameters
    ----------
    imm_seq : `<cdata 'struct imm_seq *'>`.
        Sequence pointer.
    alphabet : `T`
        Alphabet.
    """

    def __init__(self, imm_seq: ffi.CData, alphabet: T):
        if imm_seq == ffi.NULL:
            raise RuntimeError("`imm_seq` is NULL.")
        self._imm_seq = imm_seq
        self._alphabet = alphabet

    @property
    def imm_seq(self) -> ffi.CData:
        return self._imm_seq

    @property
    def length(self) -> int:
        return lib.imm_seq_length(self._imm_seq)

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_seq_string(self._imm_seq))

    def slice(self, interval: Interval) -> SubSequence[T]:
        return SubSequence[T](self, interval)

    @property
    def alphabet(self) -> T:
        return self._alphabet

    def __del__(self):
        if self._imm_seq != ffi.NULL:
            lib.imm_seq_destroy(self._imm_seq)

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Sequence(CSequence[T]):
    """
    Sequence of symbols from a given alphabet.

    Parameters
    ----------
    sequence : `bytes`
        Sequence of symbols.
    alphabet : `T`
        Alphabet.
    """

    def __init__(self, sequence: bytes, alphabet: T):
        super().__init__(lib.imm_seq_create(sequence, alphabet.imm_abc), alphabet)
        self._alphabet = alphabet

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class CSubSequence(SequenceABC[T]):
    """
    Wrapper around the C implementation of subsequence.

    Parameters
    ----------
    imm_subseq : `<cdata 'struct imm_subseq *'>`.
        Subsequence pointer.
    sequence : `CSequence`
        Sequence.
    """

    def __init__(self, imm_subseq: ffi.CData, sequence: CSequence[T]):
        if ffi.getctype(ffi.typeof(imm_subseq)) != "struct imm_subseq":
            raise TypeError("Wrong `imm_subseq` type.")
        self._imm_subseq = imm_subseq
        self._sequence = sequence

    @property
    def imm_seq(self) -> ffi.CData:
        return lib.imm_subseq_cast(ffi.addressof(self._imm_subseq))

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
        imm_seq = self.imm_seq
        return ffi.string(lib.imm_seq_string(imm_seq), lib.imm_seq_length(imm_seq))

    def slice(self, interval: Interval) -> SubSequence[T]:
        start = interval.start + self.start
        length = interval.stop - interval.start
        return SubSequence[T](self._sequence, Interval(start, start + length))

    @property
    def alphabet(self) -> T:
        return self._sequence.alphabet

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SubSequence(CSubSequence[T]):
    """
    Subsequence of symbols of a given sequence.

    Parameters
    ----------
    sequence : `CSequence`
        Sequence.
    interval : `Interval`
        Interval.
    """

    def __init__(self, sequence: CSequence[T], interval: Interval):
        length = interval.stop - interval.start
        if interval.start < 0 or length < 0 or length > sequence.length:
            raise ValueError("Out-of-range interval.")

        imm_subseq = lib.imm_subseq_slice(sequence.imm_seq, interval.start, length)
        super().__init__(imm_subseq, sequence)
