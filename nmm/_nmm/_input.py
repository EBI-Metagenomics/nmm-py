from __future__ import annotations

from typing import Type

from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import DP, HMM
from .._imm import Alphabet
from ._model import Model


class Input:
    def __init__(self, nmm_input: CData):
        if nmm_input == ffi.NULL:
            raise RuntimeError("`nmm_input` is NULL.")
        self._nmm_input = nmm_input

    @classmethod
    def create(cls: Type[Input], filepath: bytes) -> Input:
        return cls(lib.nmm_input_create(filepath))

    def read(self):
        nmm_model = lib.nmm_input_read(self._nmm_input)
        if nmm_model == ffi.NULL:
            raise RuntimeError("Could not read model.")
        abc = Alphabet(lib.nmm_model_abc(nmm_model))
        hmm = HMM(lib.nmm_model_hmm(nmm_model), abc)
        dp = DP(lib.nmm_model_dp(nmm_model), hmm)
        return Model(nmm_model, hmm, dp)

    def __del__(self):
        if self._nmm_input != ffi.NULL:
            lib.nmm_input_destroy(self._nmm_input)
