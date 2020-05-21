from __future__ import annotations

from typing import Type

from imm import DP, HMM, Alphabet

from ._cdata import CData
from ._ffi import ffi, lib

__all__ = ["Model"]


class Model:
    def __init__(self, nmm_model: CData, hmm: HMM, dp: DP):
        if nmm_model == ffi.NULL:
            raise RuntimeError("`nmm_model` is NULL.")
        self._nmm_model = nmm_model
        self._hmm = hmm
        self._dp = dp

    @property
    def nmm_model(self) -> CData:
        return self._nmm_model

    @classmethod
    def create(cls: Type[Model], hmm: HMM, dp: DP) -> Model:
        return cls(lib.nmm_model_create(hmm.imm_hmm, dp.imm_dp), hmm, dp)

    @property
    def alphabet(self) -> Alphabet:
        return self._hmm.alphabet

    @property
    def dp(self) -> DP:
        return self._dp

    @property
    def hmm(self) -> HMM:
        return self._hmm

    # @property
    # def alphabet(self) -> BaseAlphabet:
    #     return self._alphabet

    # @property
    # def nmm_codon_table(self) -> CData:
    #     return self._nmm_codon_table

    # def lprob(self, codon: Codon) -> float:
    #     lprob: float = lib.nmm_codon_table_lprob(self._nmm_codon_table, codon.nmm_codon)
    #     if not lprob_is_valid(lprob):
    #         raise RuntimeError("Could not get probability.")
    #     return lprob

    def __del__(self):
        if self._nmm_model != ffi.NULL:
            lib.nmm_model_destroy(self._nmm_model)
