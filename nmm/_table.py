from __future__ import annotations

from typing import Iterable, Tuple, Type

from imm import lprob_is_valid

from ._alphabet import AminoAlphabet, BaseAlphabet
from ._cdata import CData
from ._codon import Codon
from ._codon_prob import CodonProb
from ._ffi import ffi, lib

__all__ = ["AminoTable", "BaseTable", "CodonTable"]


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


class BaseTable:
    """
    Base table of probabilities.

    Parameters
    ----------
    nmm_base_table
        Base table.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_base_table: CData, alphabet: BaseAlphabet):
        if nmm_base_table == ffi.NULL:
            raise RuntimeError("`nmm_base_table` is NULL.")
        self._nmm_base_table = nmm_base_table
        self._alphabet = alphabet

    @classmethod
    def create(
        cls: Type[BaseTable],
        alphabet: BaseAlphabet,
        lprobs: Tuple[float, float, float, float],
    ) -> BaseTable:
        """
        Create base table of probabilities.

        Parameters
        ----------
        alphabet
            Four-nucleotides alphabet.
        lprobs
            Log probability of each nucleotide.
        """
        nmm_base_table = lib.nmm_base_table_create(alphabet.nmm_base_abc, *lprobs)
        return cls(nmm_base_table, alphabet)

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def nmm_base_table(self) -> CData:
        return self._nmm_base_table

    def lprob(self, nucleotide: bytes) -> float:
        return lib.nmm_base_table_lprob(self._nmm_base_table, nucleotide)

    def __del__(self):
        if self._nmm_base_table != ffi.NULL:
            lib.nmm_base_table_destroy(self._nmm_base_table)


class CodonTable:
    """
    Codon table.

    Compute marginal and non-marginal codon probabilities.

    Parameters
    ----------
    nmm_codon_table
        Codon table.
    alphabet
        Four-nucleotides alphabet.
    """

    def __init__(self, nmm_codon_table: CData, alphabet: BaseAlphabet):
        if nmm_codon_table == ffi.NULL:
            raise RuntimeError("`nmm_codon_table` is NULL.")
        self._nmm_codon_table = nmm_codon_table
        self._alphabet = alphabet

    @classmethod
    def create(cls: Type[CodonTable], codonp: CodonProb) -> CodonTable:
        """
        Create a codon table.

        Parameters
        ----------
        codonp
            Non-marginal codon probabilities.
        """
        return cls(lib.nmm_codon_table_create(codonp.nmm_codon_lprob), codonp.alphabet)

    @property
    def alphabet(self) -> BaseAlphabet:
        return self._alphabet

    @property
    def nmm_codon_table(self) -> CData:
        return self._nmm_codon_table

    def lprob(self, codon: Codon) -> float:
        lprob: float = lib.nmm_codon_table_lprob(self._nmm_codon_table, codon.nmm_codon)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._nmm_codon_table != ffi.NULL:
            lib.nmm_codon_table_destroy(self._nmm_codon_table)
