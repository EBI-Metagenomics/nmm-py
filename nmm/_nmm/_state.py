from typing import Dict, Tuple

from .._ffi import ffi, lib
from .._imm import Alphabet, Sequence, SequenceTable, State, TableState
from ._base_table import BaseTable
from ._codon import Codon
from ._codon_table import CodonTable
from ._base_alphabet import BaseAlphabet


class FrameState(State[BaseAlphabet]):
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
            name, baset.nmm_base_table, codont.nmm_codon_table, epsilon
        )
        if state == ffi.NULL:
            raise RuntimeError("`nmm_frame_state_create` failed.")

        self._baset = baset
        self._codont = codont
        self._epsilon = epsilon
        self._nmm_frame_state = state
        alphabet = baset.alphabet
        super().__init__(lib.imm_state_cast_c(self._nmm_frame_state), alphabet)

    def decode(self, seq: Sequence) -> Tuple[float, Codon]:
        state = self._nmm_frame_state
        any_symbol = self.alphabet.any_symbol
        codon = Codon.create(any_symbol * 3, self.alphabet)
        lprob = lib.nmm_frame_state_decode(state, seq.imm_seq, codon.nmm_codon)
        return lprob, codon

    def __del__(self):
        if self._nmm_frame_state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._nmm_frame_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class CodonState(TableState):
    def __init__(self, name: bytes, alphabet: Alphabet, emission: Dict[Codon, float]):
        """
        Parameters
        ----------
        name
            State name.
        alphabet
            Alphabet.
        emission
            Codon probabilities.
        """

        seqt = SequenceTable.create(alphabet)
        for codon, lprob in emission.items():
            seqt.add(Sequence.create(codon.symbols, codon.alphabet), lprob)

        super().__init__(name, seqt)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"
