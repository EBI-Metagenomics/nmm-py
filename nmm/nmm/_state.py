from typing import Dict

from .._ffi import ffi, lib
from ..imm import CAlphabet, CSequence, CState, Sequence, SequenceTable, TableState
from ._base_table import BaseTable
from ._codon import Codon
from ._codon_table import CodonTable


class FrameState(CState):
    def __init__(
        self, name: bytes, baset: BaseTable, codont: CodonTable, epsilon: float
    ):
        """
        Parameters
        ----------
        name : bytes
            State name.
        baset : `BaseTable`
            Base table of probabilities.
        codont : `CodonTable`
            Codon table of probabilities.
        epsilon : float
            Epsilon.
        """
        state = lib.nmm_frame_state_create(
            name, baset.nmm_baset, codont.nmm_codont, epsilon
        )
        if state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state_create` failed.")

        self._baset = baset
        self._codont = codont
        self._epsilon = epsilon
        self._nmm_frame_state = state
        alphabet = baset.base.alphabet
        super().__init__(lib.imm_state_cast_c(self._nmm_frame_state), alphabet)

    def decode(self, seq: CSequence, codon: Codon) -> float:
        state = self._nmm_frame_state
        return lib.nmm_frame_state_decode(state, seq.imm_seq, codon.nmm_codon)

    def __del__(self):
        if self._nmm_frame_state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._nmm_frame_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class CodonState(TableState):
    def __init__(self, name: bytes, alphabet: CAlphabet, emission: Dict[Codon, float]):
        """
        Parameters
        ----------
        name : bytes
            State name.
        alphabet : `CAlphabet`
            Alphabet.
        emission : `Dict[Codon, float]`
            Codon probabilities.
        """

        seqt = SequenceTable(alphabet)
        for k, v in emission.items():
            seqt.add(Sequence(k.symbols, k.base.alphabet), v)

        super().__init__(name, seqt)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
