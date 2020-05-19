from enum import Enum
from .._cdata import CData
from ._base_alphabet import BaseAlphabet
from ._amino_alphabet import AminoAlphabet
from .._imm import Alphabet

from .._ffi import lib


class AlphabetType(Enum):
    BASE = 0x10
    AMINO = 0x11


def wrap_imm_abc(imm_abc: CData):
    try:
        alphabet_type = AlphabetType(lib.imm_abc_type_id(imm_abc))
    except ValueError:
        return Alphabet(imm_abc)
    if alphabet_type == AlphabetType.BASE:
        nmm_base_abc = lib.nmm_base_abc_derived(imm_abc)
        return BaseAlphabet(nmm_base_abc)
    if alphabet_type == AlphabetType.AMINO:
        nmm_amino_abc = lib.nmm_amino_abc_derived(imm_abc)
        return AminoAlphabet(nmm_amino_abc)
    raise RuntimeError("It should not get here.")
