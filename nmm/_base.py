from ._alphabet import Alphabet
from ._string import make_sure_bytes

from ._ffi import ffi, lib


class Base:
    def __init__(self, alphabet: Alphabet, lprobs: dict = {}):
        self._alphabet = alphabet
        self._base = ffi.NULL
        self._base = lib.nmm_base_create(self._alphabet.cdata)
        for letter, lprob in lprobs.items():
            self.set_lprob(letter, lprob)

    @property
    def cdata(self):
        return self._base

    @property
    def alphabet(self):
        return self._alphabet

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
