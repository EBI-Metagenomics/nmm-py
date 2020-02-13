from ._lprob import lprob_is_valid
from ._alphabet import CAlphabet
from ._sequence import CSequence
from .._ffi import ffi, lib


class CSequenceTable:
    """
    Wrapper around the C implementation of sequence table.

    Parameters
    ----------
    imm_seq_table : `<cdata 'struct imm_seq_table *'>`.
    """

    def __init__(self, imm_seq_table: ffi.CData):
        if imm_seq_table == ffi.NULL:
            raise RuntimeError("`imm_seq_table` is NULL.")
        self._imm_seq_table = imm_seq_table

    @property
    def imm_seq_table(self) -> ffi.CData:
        return self._imm_seq_table

    @property
    def length(self) -> int:
        imm_abc = lib.imm_seq_table_get_abc(self._imm_seq_table)
        return lib.imm_abc_length(imm_abc)

    @property
    def symbols(self) -> bytes:
        imm_abc = lib.imm_seq_table_get_abc(self._imm_seq_table)
        return ffi.string(lib.imm_abc_string(imm_abc))

    def add(self, sequence: CSequence, lprob: float):
        if lib.imm_seq_table_add(self._imm_seq_table, sequence.imm_seq, lprob) != 0:
            raise RuntimeError("Could not add sequence.")

    def normalize(self):
        if lib.imm_seq_table_normalize(self._imm_seq_table) != 0:
            raise RuntimeError("Could not normalize it.")

    def lprob(self, sequence: CSequence) -> float:
        lprob: float = lib.imm_seq_table_lprob(self._imm_seq_table, sequence.imm_seq)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __del__(self):
        if self._imm_seq_table != ffi.NULL:
            lib.imm_seq_table_destroy(self._imm_seq_table)

    def __str__(self) -> str:
        return f"[{self.symbols.decode()}]"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class SequenceTable(CSequenceTable):
    """
    Table of sequence probabilities.

    Parameters
    ----------
    alphabet : `CAlphabet`
        Alphabet.
    """

    def __init__(self, alphabet: CAlphabet):
        self._alphabet = alphabet
        super().__init__(lib.imm_seq_table_create(alphabet.imm_abc))

    @property
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
