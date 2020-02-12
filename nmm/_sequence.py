from abc import ABC, abstractmethod

from ._alphabet import CAlphabet
from ._ffi import ffi, lib


class SequenceABC(ABC):
    @property
    @abstractmethod
    def length(self) -> int:
        raise NotImplementedError()

    @property
    @abstractmethod
    def symbols(self) -> bytes:
        raise NotImplementedError()


class CSequence(SequenceABC):
    """
    Wrapper around the C implementation of sequence.

    Parameters
    ----------
    imm_seq : `<cdata 'struct imm_seq *'>`.
        Sequence pointer.
    """

    def __init__(self, imm_seq: ffi.CData):
        if imm_seq == ffi.NULL:
            raise RuntimeError("`imm_seq` is NULL.")
        self._imm_seq = imm_seq

    @property
    def imm_seq(self) -> ffi.CData:
        return self._imm_seq

    @property
    def length(self) -> int:
        return lib.imm_seq_length(self._imm_seq)

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_seq_string(self._imm_seq))

    def __del__(self):
        if self._imm_seq != ffi.NULL:
            lib.imm_seq_destroy(self._imm_seq)

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Sequence(CSequence):
    """
    Sequence of symbols from a given alphabet.

    Parameters
    ----------
    seq : bytes
        Sequence of symbols.
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, seq: bytes, alphabet: CAlphabet):
        self._alphabet = alphabet
        super().__init__(lib.imm_seq_create(seq, alphabet.imm_abc))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
