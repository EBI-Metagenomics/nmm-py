from abc import ABC, abstractmethod
from typing import Tuple

from ._base import Base
from ._ffi import ffi, lib


class CodonABC(ABC):
    @abstractmethod
    def get(self) -> Tuple[int, int, int]:
        raise NotImplementedError()

    @abstractmethod
    def set(self, codon: Tuple[int, int, int]):
        del codon
        raise NotImplementedError()

    # @property
    # @abstractmethod
    # def nmm_base(self) -> ffi.CData:
    #     raise NotImplementedError()

    @property
    @abstractmethod
    def nmm_codon(self) -> ffi.CData:
        raise NotImplementedError()


class CCodon(CodonABC):
    def __init__(self, nmm_codon: ffi.CData):
        super().__init__()
        if nmm_codon == ffi.NULL:
            raise RuntimeError("`nmm_codon` is NULL.")
        self._nmm_codon = nmm_codon

    def get(self) -> Tuple[int, int, int]:
        triplet = lib.nmm_codon_get(self._nmm_codon)
        return (triplet.a, triplet.b, triplet.c)

    def set(self, codon: Tuple[int, int, int]):
        triplet = ffi.new("struct nmm_triplet *")
        triplet.a = codon[0]
        triplet.b = codon[1]
        triplet.c = codon[2]
        e: int = lib.nmm_codon_set(self._nmm_codon, triplet)
        if e != 0:
            raise ValueError("Could not set codon")

    # @property
    # def nmm_base(self) -> ffi.CData:
    #     return lib.nmm_codon_get_base(self._nmm_codon)

    @property
    def nmm_codon(self) -> ffi.CData:
        return self._nmm_codon

    def __del__(self):
        if self._nmm_codon != ffi.NULL:
            lib.nmm_codon_destroy(self._nmm_codon)


class Codon(CCodon):
    """
    Codon is a sequence of three bases.
    """

    def __init__(self, codon: Tuple[int, int, int], base: Base):
        super().__init__(lib.nmm_codon_create(base.nmm_base))
        self.set(codon)
        self._base = base

    def __eq__(self, another):
        return bytes(self) == bytes(another)

    def __hash__(self):
        return hash(bytes(self))

    def __bytes__(self) -> bytes:
        triplet = self.get()
        a = chr(triplet[0]).encode()
        b = chr(triplet[1]).encode()
        c = chr(triplet[2]).encode()
        return a + b + c

    def __str__(self) -> str:
        return bytes(self).decode()

    def __repr__(self) -> str:
        codon = bytes(self)
        return f"<{self.__class__.__name__}:{codon.decode()}>"
