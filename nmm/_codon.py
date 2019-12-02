from ._alphabet import Alphabet, CAlphabet
from typing import Dict, Optional, Union
from ._base import Base
from ._ffi import ffi, lib


class Codon:
    """
    Codon is a sequence of three bases.

    codon : Union[bytes, str, ffi.CData]
        Sequence of letters representing a codon.
    """

    def __init__(self, codon: Union[bytes, str, ffi.CData]):
        if isinstance(codon, str):
            codon = codon.encode()

        if isinstance(codon, bytes):
            if len(codon) != 3:
                raise ValueError("Codon must have three letters.")

            nmm_codon = ffi.new("struct nmm_codon *")
            nmm_codon.a = bytes(codon)[0:1]
            nmm_codon.b = bytes(codon)[1:2]
            nmm_codon.c = bytes(codon)[2:3]
            self._nmm_codon = nmm_codon
        else:
            self._nmm_codon = codon

    @property
    def nmm_codon(self) -> ffi.CData:
        return self._nmm_codon

    def base(self, idx: int) -> Base:
        if idx == 0:
            return Base(self._nmm_codon.a)
        elif idx == 1:
            return Base(self._nmm_codon.b)
        elif idx == 2:
            return Base(self._nmm_codon.c)
        raise ValueError("Index must be 0, 1, or 2.")

    def __eq__(self, another):
        return bytes(self) == bytes(another)

    def __hash__(self):
        return hash(bytes(self))

    def __bytes__(self) -> bytes:
        return self._nmm_codon.a + self._nmm_codon.b + self._nmm_codon.c

    def __str__(self) -> str:
        return bytes(self).decode()

    def __repr__(self) -> str:
        codon = bytes(self)
        return f"<{self.__class__.__name__}:{codon.decode()}>"


class CCodonTable:
    """
    Wrapper around the C implementation of a codon table.

    Parameters
    ----------
    nmm_codont : Optional[ffi.CData]
        Passing `None` will create a new codon table at the underlying library level using the
        `alphabet` argument.
    alphabet : Optional[CAlphabet]
        Passing a `CAlphabet` will create a new codon table at the underlying library level.
    """

    def __init__(
        self,
        nmm_codont: Optional[ffi.CData] = None,
        alphabet: Optional[CAlphabet] = None,
    ):
        self._nmm_codont = ffi.NULL
        if nmm_codont is None:
            if alphabet is None:
                raise ValueError("`alphabet` is `None`")

            self._nmm_codont = lib.nmm_codont_create(alphabet.imm_abc)
        else:
            if alphabet is not None:
                raise ValueError("`alphabet` is not `None`")

            self._nmm_codont = nmm_codont

        if self._nmm_codont == ffi.NULL:
            raise RuntimeError("`cdata` is NULL.")

    @property
    def nmm_codont(self) -> ffi.CData:
        return self._nmm_codont

    @property
    def imm_abc(self) -> ffi.CData:
        return lib.nmm_codont_get_abc(self._nmm_codont)

    def set_lprob(self, codon: Codon, lprob: float) -> None:
        err: int = lib.nmm_codont_set_lprob(self._nmm_codont, codon.nmm_codon, lprob)
        if err != 0:
            raise ValueError(f"Could not set a probability for `{str(codon)}`.")

    def get_lprob(self, codon: Codon) -> float:
        return lib.nmm_codont_get_lprob(self._nmm_codont, codon.nmm_codon)

    def normalize(self) -> None:
        err: int = lib.nmm_codont_normalize(self._nmm_codont)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __del__(self):
        if self._nmm_codont != ffi.NULL:
            lib.nmm_codont_destroy(self._nmm_codont)


class CodonTable(CCodonTable):
    """
    Codon table.

    It is a table of probabilities for codons, in log-space.

    Parameters
    ----------
    alphabet : Alphabet.
    lprobs : Emission probabilities in log-space.
    """

    def __init__(self, alphabet: Alphabet, lprobs: Dict[Codon, float] = {}):
        super().__init__(alphabet=alphabet)
        self._alphabet = alphabet
        for seq, lprob in lprobs.items():
            self.set_lprob(seq, lprob)

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet
