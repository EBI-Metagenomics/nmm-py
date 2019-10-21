from ._alphabet import Alphabet
from ._string import make_sure_bytes

from ._ffi import ffi, lib


class Base:
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._base = ffi.NULL
        self._base = lib.nmm_base_create(self._alphabet.cdata)

    def set_lprob(self, nucleotide: str, lprob: float):
        nucleotide = make_sure_bytes(nucleotide)
        if len(nucleotide) != 1:
            raise ValueError("Nucleotide must be a single letter.")
        err: int = lib.nmm_base_set_lprob(self._base, nucleotide, lprob)
        if err != 0:
            raise ValueError(f"Could not set a probability for `{nucleotide}`.")

    def get_lprob(self, nucleotide: str) -> float:
        nucleotide = make_sure_bytes(nucleotide)
        if len(nucleotide) != 1:
            raise ValueError("Nucleotide must be a single letter.")
        return lib.nmm_base_get_lprob(self._base, nucleotide)

    def normalize(self):
        err: int = lib.nmm_base_normalize(self._base)
        if err != 0:
            raise ValueError("Normalization error.")

    def __del__(self):
        if self._base != ffi.NULL:
            lib.nmm_base_destroy(self._base)
