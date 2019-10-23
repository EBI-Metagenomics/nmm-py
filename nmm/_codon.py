from ._alphabet import Alphabet
from typing import Dict

from ._ffi import ffi, lib


class Codon:
    def __init__(self, alphabet: Alphabet, lprobs: Dict[str, float] = {}):
        self._alphabet = alphabet
        self._codon = ffi.NULL
        self._codon = lib.nmm_codon_create(self._alphabet.cdata)
        for seq, lprob in lprobs.items():
            self.set_lprob(seq, lprob)

    @property
    def cdata(self) -> ffi.CData:
        return self._codon

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    def set_lprob(self, seq: str, lprob: float) -> None:
        s = seq.encode()
        if len(s) != 3:
            raise ValueError("Codon must have three letters.")
        err: int = lib.nmm_codon_set_lprob(self._codon, s[0:1], s[1:2], s[2:3], lprob)
        if err != 0:
            raise ValueError(f"Could not set a probability for `{seq}`.")

    def get_lprob(self, seq: str) -> float:
        s = seq.encode()
        if len(s) != 3:
            raise ValueError("Codon must have three letters.")
        return lib.nmm_codon_get_lprob(self._codon, s[0:1], s[1:2], s[2:3])

    def normalize(self) -> None:
        err = lib.nmm_codon_normalize(self._codon)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self._codon != ffi.NULL:
            lib.nmm_codon_destroy(self._codon)
