from ._alphabet import Alphabet
from ._string import make_sure_bytes

from ._ffi import ffi, lib


class Codon:
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._codon = ffi.NULL
        self._codon = lib.nmm_codon_create(self._alphabet.cdata)

    def set_lprob(self, seq: str, lprob: float):
        seq = make_sure_bytes(seq)
        lib.nmm_codon_set_lprob(self._codon, seq[0:1], seq[1:2], seq[2:3], lprob)

    def get_lprob(self, seq: str) -> float:
        seq = make_sure_bytes(seq)
        return lib.nmm_codon_get_lprob(self._codon, seq[0:1], seq[1:2], seq[2:3])

    def normalize(self):
        err = lib.nmm_codon_normalize(self._codon)
        if err != 0:
            raise ValueError("Normalization error.")

    def __del__(self):
        if self._codon != ffi.NULL:
            lib.nmm_codon_destroy(self._codon)
