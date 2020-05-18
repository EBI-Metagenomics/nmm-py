from __future__ import annotations

from typing import Mapping, Type

from .._cdata import CData
from .._ffi import ffi, lib
from .._imm import DP, HMM, Alphabet
from ._model import Model
from ._state import State, wrap_imm_state


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
        nstates: int = lib.nmm_model_nstates(nmm_model)
        states: Mapping[CData, State] = {}
        for i in range(nstates):
            imm_state = lib.nmm_model_state(nmm_model, i)
            states[imm_state] = wrap_imm_state(imm_state, abc)

        hmm = HMM(lib.nmm_model_hmm(nmm_model), abc, states)
        dp = DP(lib.nmm_model_dp(nmm_model), hmm)
        return Model(nmm_model, hmm, dp)

    def close(self):
        err: int = lib.nmm_input_close(self._nmm_input)
        if err != 0:
            raise RuntimeError("Could not close input.")

    def __del__(self):
        if self._nmm_input != ffi.NULL:
            lib.nmm_input_destroy(self._nmm_input)
