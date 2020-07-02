from __future__ import annotations

from typing import Dict, Iterator, Type

from imm import DP, HMM, State

from . import wrap
from ._alphabet import BaseAlphabet
from ._cdata import CData
from ._codon_prob import CodonProb
from ._ffi import ffi, lib
from ._model import Model
from ._table import BaseTable, CodonTable

__all__ = ["Input"]


class Input:
    def __init__(self, nmm_input: CData):
        if nmm_input == ffi.NULL:
            raise RuntimeError("`nmm_input` is NULL.")
        self._nmm_input = nmm_input

    @classmethod
    def create(cls: Type[Input], filepath: bytes) -> Input:
        return cls(lib.nmm_input_create(filepath))

    def fseek(self, offset: int):
        err: int = lib.nmm_input_fseek(self._nmm_input, offset)
        if err != 0:
            raise RuntimeError("Could not fseek.")

    def ftell(self) -> int:
        offset: int = lib.nmm_input_ftell(self._nmm_input)
        if offset < 0:
            raise RuntimeError("Could not ftell.")
        return offset

    def read(self) -> Model:
        nmm_model = lib.nmm_input_read(self._nmm_input)
        if nmm_model == ffi.NULL:
            if lib.nmm_input_eof(self._nmm_input):
                raise StopIteration
            raise RuntimeError("Could not read model.")

        abc = wrap.imm_abc(lib.nmm_model_abc(nmm_model))

        base_tables = read_base_tables(nmm_model, abc)
        codon_tables = read_codon_tables(nmm_model, abc)
        codon_probs = read_codon_probs(nmm_model, abc)

        states: Dict[CData, State] = {}
        for i in range(lib.nmm_model_nstates(nmm_model)):
            imm_state = lib.nmm_model_state(nmm_model, i)
            states[imm_state] = wrap.imm_state(
                imm_state, abc, base_tables, codon_tables, codon_probs
            )

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

    def __iter__(self) -> Iterator[Model]:
        while True:
            try:
                yield self.read()
            except StopIteration:
                return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()


def read_base_tables(nmm_model: CData, abc: BaseAlphabet) -> Dict[CData, BaseTable]:
    base_tables: Dict[CData, BaseTable] = {}
    for i in range(lib.nmm_model_nbase_tables(nmm_model)):
        nmm_base_table = lib.nmm_model_base_table(nmm_model, i)
        base_tables[nmm_base_table] = BaseTable(nmm_base_table, abc)
    return base_tables


def read_codon_tables(nmm_model: CData, abc: BaseAlphabet) -> Dict[CData, CodonTable]:
    codon_tables: Dict[CData, CodonTable] = {}
    for i in range(lib.nmm_model_ncodon_tables(nmm_model)):
        nmm_codon_table = lib.nmm_model_codon_table(nmm_model, i)
        codon_tables[nmm_codon_table] = CodonTable(nmm_codon_table, abc)
    return codon_tables


def read_codon_probs(nmm_model: CData, abc: BaseAlphabet) -> Dict[CData, CodonProb]:
    codon_probs: Dict[CData, CodonProb] = {}
    for i in range(lib.nmm_model_ncodon_lprobs(nmm_model)):
        nmm_codon_lprob = lib.nmm_model_codon_lprob(nmm_model, i)
        codon_probs[nmm_codon_lprob] = CodonProb(nmm_codon_lprob, abc)
    return codon_probs
