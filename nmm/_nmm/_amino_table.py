from typing import Iterable

from .._ffi import ffi, lib
from ._amino_alphabet import CAminoAlphabet


class CAminoTable:
    """
    Wrapper around the C implementation of a amino table.

    Parameters
    ----------
    nmm_amino_table : `<cdata 'struct nmm_amino_table *'>`.
        Amino table.
    amino : `CAminoAlphabet`
        20-nucleotides alphabet.
    """

    def __init__(self, nmm_amino_table: ffi.CData, amino_abc: CAminoAlphabet):
        if nmm_amino_table == ffi.NULL:
            raise RuntimeError("`nmm_amino_table` is NULL.")
        self._nmm_amino_table = nmm_amino_table
        self._amino_abc = amino_abc

    @property
    def alphabet(self) -> CAminoAlphabet:
        return self._amino_abc

    @property
    def nmm_amino_table(self) -> ffi.CData:
        return self._nmm_amino_table

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_amino_table_lprob(self._nmm_amino_table, nucleotide)

    def __del__(self):
        if self._nmm_amino_table != ffi.NULL:
            lib.nmm_amino_table_destroy(self._nmm_amino_table)


class AminoTable(CAminoTable):
    """
    Amino table of probabilities.

    Parameters
    ----------
    amino : `CAminoAlphabet`
        20-nucleotides alphabet.
    lprobs : `Tuple[float, float, float, float]`
        Log probability of each nucleotide.
    """

    def __init__(self, amino_abc: CAminoAlphabet, lprobs: Iterable[float]):
        nmm_amino_table = lib.nmm_amino_table_create(
            amino_abc.nmm_amino_abc, list(lprobs)
        )
        super().__init__(nmm_amino_table, amino_abc)
