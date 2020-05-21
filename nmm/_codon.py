from __future__ import annotations

import itertools
from typing import Iterable, Type

from ._alphabet import BaseAlphabet
from ._cdata import CData
from ._ffi import ffi, lib

__all__ = [
    "Codon",
    "codon_iter",
]


class Codon:
    """
    Codon is a sequence of three symbols from a four-nucleotides alphabet.

    Parameters
    ----------
    nmm_codon
        Codon pointer.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon: CData, alphabet: BaseAlphabet):
        if nmm_codon == ffi.NULL:
            raise RuntimeError("`nmm_codon` is NULL.")
        self._nmm_codon = nmm_codon
        self._alphabet = alphabet

    @classmethod
    def create(cls: Type[Codon], symbols: bytes, alphabet: BaseAlphabet) -> Codon:
        """
        Create a codon.

        Parameters
        ----------
        symbols
            Sequence of four symbols.
        alphabet
            Four-nucleotides alphabet.
        """
        codon = cls(lib.nmm_codon_create(alphabet.nmm_base_abc), alphabet)
        codon.symbols = symbols
        return codon

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def symbols(self) -> bytes:
        triplet = lib.nmm_codon_get_triplet(self._nmm_codon)
        return triplet.a + triplet.b + triplet.c

    @symbols.setter
    def symbols(self, symbols: bytes):
        if len(symbols) != 3:
            raise ValueError("Symbols length must be three.")

        triplet = {"a": symbols[0:1], "b": symbols[1:2], "c": symbols[2:3]}
        if lib.nmm_codon_set_triplet(self._nmm_codon, triplet) != 0:
            raise ValueError("Could not set codon.")

    @property
    def nmm_codon(self) -> CData:
        return self._nmm_codon

    def __del__(self):
        if self._nmm_codon != ffi.NULL:
            lib.nmm_codon_destroy(self._nmm_codon)

    def __eq__(self, another):
        return bytes(self) == bytes(another)

    def __hash__(self):
        return hash(bytes(self))

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __bytes__(self) -> bytes:
        return str(self).encode()

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


def codon_iter(base_abc: BaseAlphabet) -> Iterable[Codon]:
    """
    Codon iterator.

    Parameters
    ----------
    base_abc
        Base alphabet.
    """
    bases = [base_abc.symbols[i : i + 1] for i in range(len(base_abc.symbols))]

    for a, b, c in itertools.product(bases, bases, bases):
        yield Codon.create(a + b + c, base_abc)
