from itertools import product
from scipy.special import logsumexp
from math import exp

from ._string import make_sure_bytes
from ._log import LOG
from ._norm import normalize_emission
from ._codon import Codon
from ._state import State, emission_table

from ._ffi import ffi, lib

# struct nmm_frame_state *nmm_frame_state_create(const char *name, const struct imm_abc *bases,
#                                                const double *          base_lprobs,
#                                                const struct nmm_codon *codon, double epsilon);
# void                    nmm_frame_state_destroy(struct nmm_frame_state *state);
# int                     nmm_frame_state_normalize(struct nmm_frame_state *state);


class FrameState(State):
    def __init__(self, name: str, codon: Codon, emission: dict):
        """
        Parameters
        ----------
        name : str
            Name.
        alphabet : Alphabet
            Alphabet.
        emission : dict
            Emission probabilities in log-space.
        """
        super(FrameState, self).__init__(codon.alphabet)

        name = make_sure_bytes(name)
        self._state = lib.imm_table_state_create(name, codon.alphabet.cdata)

        for seq, lprob in emission.items():
            seq = make_sure_bytes(seq)
            lib.imm_table_state_add(self._state, seq, lprob)

    def normalize(self):
        err = lib.imm_table_state_normalize(self._state)
        if err != 0:
            raise ValueError("Normalization error.")

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_table_state_destroy(self._state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"
