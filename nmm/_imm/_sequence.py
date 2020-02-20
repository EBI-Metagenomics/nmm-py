from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Generic, Type, TypeVar, Union

from .._ffi import ffi, lib
from .._cdata import CData
from .._interval import Interval
from ._alphabet import Alphabet

T = TypeVar("T", bound=Alphabet)


class SequenceABC(Generic[T], ABC):
    """
    Sequence of symbols.
    """

    @abstractmethod
    def __len__(self) -> int:
        raise NotImplementedError()

    @abstractmethod
    def __bytes__(self) -> bytes:
        raise NotImplementedError()

    @abstractmethod
    def __getitem__(self, i: Union[int, slice, Interval]):
        del i
        raise NotImplementedError()

    @property
    @abstractmethod
    def alphabet(self) -> T:
        raise NotImplementedError()


class Sequence(SequenceABC[T]):
    """
    Sequence of symbols from a given alphabet.

    Parameters
    ----------
    imm_seq
        Sequence pointer.
    alphabet
        Alphabet.
    """

    def __init__(self, imm_seq: CData, alphabet: T):
        if imm_seq == ffi.NULL:
            raise RuntimeError("`imm_seq` is NULL.")
        self._imm_seq = imm_seq
        self._alphabet = alphabet

    @classmethod
    def create(cls: Type[Sequence[T]], sequence: bytes, alphabet: T) -> Sequence[T]:
        """
        Create a sequence of symbols.

        Parameters
        ----------
        sequence
            Sequence of symbols.
        alphabet
            Alphabet.
        """
        return cls(lib.imm_seq_create(sequence, alphabet.imm_abc), alphabet)

    @property
    def imm_seq(self) -> CData:
        return self._imm_seq

    def __len__(self) -> int:
        return lib.imm_seq_length(self._imm_seq)

    def __bytes__(self) -> bytes:
        return ffi.string(lib.imm_seq_string(self._imm_seq))

    def __getitem__(self, i: Union[int, slice, Interval]):
        if isinstance(i, int):
            return bytes(self)[i : i + 1]
        elif isinstance(i, slice):
            interval = Interval.from_slice(i)
        elif isinstance(i, Interval):
            interval = i
        else:
            raise RuntimeError("Index has to be an integer of a slice.")

        return SubSequence[T].create(self, interval)

    @property
    def alphabet(self) -> T:
        return self._alphabet

    def __del__(self):
        if self._imm_seq != ffi.NULL:
            lib.imm_seq_destroy(self._imm_seq)

    def __str__(self) -> str:
        return f"{bytes(self).decode()}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SubSequence(SequenceABC[T]):
    """
    Subsequence of symbols of a given sequence.

    Parameters
    ----------
    imm_subseq
        Subsequence pointer.
    sequence
        Sequence.
    """

    def __init__(self, imm_subseq: CData, sequence: Sequence[T]):
        if ffi.getctype(ffi.typeof(imm_subseq)) != "struct imm_subseq":
            raise TypeError("Wrong `imm_subseq` type.")
        self._imm_subseq = imm_subseq
        self._sequence = sequence

    @classmethod
    def create(
        cls: Type[SubSequence[T]], sequence: Sequence[T], interval: Interval
    ) -> SubSequence[T]:
        """
        Subsequence of symbols of a given sequence.

        Parameters
        ----------
        sequence
            Sequence.
        interval
            Interval.
        """
        length = interval.stop - interval.start
        if interval.start < 0 or length < 0 or length > len(sequence):
            raise ValueError("Out-of-range interval.")

        imm_subseq = lib.imm_subseq_slice(sequence.imm_seq, interval.start, length)
        return cls(imm_subseq, sequence)

    @property
    def imm_seq(self) -> CData:
        return lib.imm_subseq_cast(ffi.addressof(self._imm_subseq))

    @property
    def imm_subseq(self) -> CData:
        return self._imm_subseq

    @property
    def start(self) -> int:
        return lib.imm_subseq_start(ffi.addressof(self._imm_subseq))

    def __len__(self) -> int:
        return lib.imm_subseq_length(ffi.addressof(self._imm_subseq))

    def __bytes__(self) -> bytes:
        imm_seq = self.imm_seq
        return ffi.string(lib.imm_seq_string(imm_seq), lib.imm_seq_length(imm_seq))

    def __getitem__(self, i: Union[int, slice, Interval]):
        if isinstance(i, int):
            return bytes(self)[i : i + 1]
        elif isinstance(i, slice):
            interval = Interval.from_slice(i)
        elif isinstance(i, Interval):
            interval = i
        else:
            raise RuntimeError("Index has to be an integer of a slice.")

        start = interval.start + self.start
        length = interval.stop - interval.start
        return SubSequence[T].create(self._sequence, Interval(start, start + length))

    @property
    def alphabet(self) -> T:
        return self._sequence.alphabet

    def __str__(self) -> str:
        return f"{bytes(self).decode()}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
