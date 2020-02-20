from __future__ import annotations

from typing import Iterable, Type

from .._cdata import CData
from .._ffi import ffi, lib
from ._amino_alphabet import AminoAlphabet


class AminoTable:
    """
    Amino table of probabilities.

    Parameters
    ----------
    nmm_amino_table
        Amino table.
    alphabet
        20-symbols alphabet.
    """

    def __init__(self, nmm_amino_table: CData, alphabet: AminoAlphabet):
        if nmm_amino_table == ffi.NULL:
            raise RuntimeError("`nmm_amino_table` is NULL.")
        self._nmm_amino_table = nmm_amino_table
        self._alphabet = alphabet

    @classmethod
    def create(
        cls: Type[AminoTable], alphabet: AminoAlphabet, lprobs: Iterable[float],
    ) -> AminoTable:
        """
        Create an amino table of probabilities.

        Parameters
        ----------
        alphabet
            20-symbols alphabet.
        lprobs
            Log probability of each amino acid.
        """
        abc = alphabet.nmm_amino_abc
        nmm_amino_table = lib.nmm_amino_table_create(abc, list(lprobs))
        return cls(nmm_amino_table, alphabet)

    @property
    def alphabet(self) -> AminoAlphabet:
        return self._alphabet

    @property
    def nmm_amino_table(self) -> CData:
        return self._nmm_amino_table

    def lprob(self, amino: bytes) -> float:
        return lib.nmm_amino_table_lprob(self._nmm_amino_table, amino)

    def __del__(self):
        if self._nmm_amino_table != ffi.NULL:
            lib.nmm_amino_table_destroy(self._nmm_amino_table)
