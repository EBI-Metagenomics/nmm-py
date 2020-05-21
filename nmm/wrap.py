from typing import Dict, TypeVar

import imm

from ._alphabet import AlphabetType, AminoAlphabet, BaseAlphabet
from ._cdata import CData
from ._codon_prob import CodonProb
from ._ffi import ffi, lib
from ._state import CodonState, FrameState, StateType
from ._table import BaseTable, CodonTable

T = TypeVar("T", bound=imm.Alphabet)


def imm_abc(ptr: CData):
    try:
        alphabet_type = AlphabetType(imm.lib.imm_abc_type_id(ptr))
    except ValueError:
        return imm.Alphabet(ptr)

    if alphabet_type == AlphabetType.BASE:
        nmm_base_abc = lib.nmm_base_abc_derived(ptr)
        return BaseAlphabet(nmm_base_abc)

    if alphabet_type == AlphabetType.AMINO:
        nmm_amino_abc = lib.nmm_amino_abc_derived(ptr)
        return AminoAlphabet(nmm_amino_abc)

    raise RuntimeError("It should not get here.")


def imm_state(
    ptr: CData,
    alphabet: T,
    base_tables: Dict[CData, BaseTable],
    codon_tables: Dict[CData, CodonTable],
    codon_probs: Dict[CData, CodonProb],
) -> imm.State:
    try:
        state_type = StateType(imm.lib.imm_state_type_id(ptr))
    except ValueError:
        return imm.wrap.imm_state(ptr, alphabet)

    if state_type == StateType.CODON:
        nmm_codon_state = lib.nmm_codon_state_derived(ptr)
        if nmm_codon_state == ffi.NULL:
            raise RuntimeError("`nmm_codon_state` is NULL.")

        nmm_codon_lprob = lib.nmm_codon_state_codon_lprob(nmm_codon_state)
        codonp = codon_probs[nmm_codon_lprob]
        return CodonState(nmm_codon_lprob, codonp)

    if state_type == StateType.FRAME:
        nmm_frame_state = lib.nmm_frame_state_derived(ptr)
        if nmm_frame_state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state` is NULL.")

        nmm_base_table = lib.nmm_frame_state_base_table(nmm_frame_state)
        nmm_codon_table = lib.nmm_frame_state_codon_table(nmm_frame_state)

        baset = base_tables[nmm_base_table]
        codont = codon_tables[nmm_codon_table]

        return FrameState(nmm_frame_state, baset, codont)

    raise ValueError(f"Unknown state type: {imm.lib.imm_state_type_id(ptr)}.")
