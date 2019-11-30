from typing import Dict, NamedTuple, Tuple

from ._alphabet import Alphabet
from ._base import BaseTable
from ._codon import CodonTable
from ._ffi import ffi, lib
from ._log import LOG0

DecodedCodon = NamedTuple("DecodedCodon", [("lprob", float), ("codon", bytes)])


# TODO: consider removing CState (and other C classes) as base classes,
# and instead to help compose their respective Python classes.
class CState:
    def __init__(self, cdata: ffi.CData):
        self.__cdata = cdata

    @property
    def name(self) -> bytes:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self.__cdata == ffi.NULL:
            raise RuntimeError("State has failed to initialize.")
        return ffi.string(lib.imm_state_get_name(self.__cdata))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self.__cdata)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self.__cdata)

    @property
    def imm_state(self) -> ffi.CData:
        return self.__cdata

    def lprob(self, seq: bytes) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        seq : Sequence.
        """
        return lib.imm_state_lprob(self.__cdata, seq, len(seq))


class State(CState):
    def __init__(self, cdata: ffi.CData, alphabet: Alphabet):
        """
        Parameters
        ----------
        cdata : `struct imm_state *`
        alphabet : Alphabet.
        """
        super().__init__(cdata)
        self._alphabet = alphabet

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    def __str__(self) -> str:
        return f"<{self.name.decode()}>"


class MuteState(State):
    def __init__(self, name: bytes, alphabet: Alphabet):
        """
        Parameters
        ----------
        name : Name.
        alphabet : Alphabet.
        """
        cdata = lib.imm_mute_state_create(name, alphabet.imm_abc)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")
        self._cdata = cdata
        super().__init__(lib.imm_state_cast_c(cdata), alphabet)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.imm_mute_state_destroy(self._cdata)


class NormalState(State):
    def __init__(self, name: bytes, alphabet: Alphabet, lprobs: Dict[bytes, float]):
        """
        Parameters
        ----------
        name : Name.
        alphabet : Alphabet.
        lprobs : Emission probabilities in log-space.
        """

        if len(set(b"".join(lprobs.keys())) - set(alphabet.symbols)) > 0:
            raise ValueError("Unrecognized alphabet symbol.")

        arr = [lprobs.get(bytes([symb]), LOG0) for symb in alphabet.symbols]
        cdata = lib.imm_normal_state_create(name, alphabet.imm_abc, arr)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        self._cdata = cdata
        super().__init__(lib.imm_state_cast_c(cdata), alphabet)

    def emission_table(self) -> Dict[bytes, float]:
        return {bytes([s]): self.lprob(bytes([s])) for s in self.alphabet.symbols}

    def normalize(self) -> None:
        err = lib.imm_normal_state_normalize(self._cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.imm_normal_state_destroy(self._cdata)


class TableState(State):
    def __init__(self, name: bytes, alphabet: Alphabet, emission: Dict[bytes, float]):
        """
        Parameters
        ----------
        name : Name.
        alphabet : Alphabet.
        emission : Emission probabilities in log-space.
        """
        cdata = lib.imm_table_state_create(name, alphabet.imm_abc)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_table_state_create` failed.")

        for seq, lprob in emission.items():
            lib.imm_table_state_add(cdata, seq, lprob)

        self._cdata = cdata
        super().__init__(lib.imm_state_cast_c(cdata), alphabet)

    def normalize(self) -> None:
        err = lib.imm_table_state_normalize(self._cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.imm_table_state_destroy(self._cdata)


class CodonState(TableState):
    # TODO: consider creating a Codon type and place it in the keys as a type.
    def __init__(self, name: bytes, alphabet: Alphabet, emission: Dict[bytes, float]):
        if sum(len(k) != 3 for k in emission.keys()) > 0:
            raise ValueError("Codon must be composed of three bases.")

        super().__init__(name, alphabet, emission)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"


class FrameState(State):
    def __init__(self, name: bytes, base: BaseTable, codon: CodonTable, epsilon: float):
        """
        Parameters
        ----------
        name : Name.
        base : Base table.
        codon : Codon table.
        epsilon : Epsilon.
        """
        if set(base.alphabet.symbols) != set(codon.alphabet.symbols):
            raise ValueError("Alphabet symbols of `base` and `codon` are not equal.")

        self._base = base
        self._codon = codon
        self._epsilon = epsilon

        cdata = lib.nmm_frame_state_create(
            name, base.nmm_baset, codon.nmm_codont, epsilon
        )
        if cdata == ffi.NULL:
            raise RuntimeError("Could not create state.")

        self._cdata = cdata
        super().__init__(lib.imm_state_cast_c(cdata), codon.alphabet)

    @property
    def base(self) -> BaseTable:
        return self._base

    @property
    def codon(self) -> CodonTable:
        return self._codon

    @property
    def epsilon(self):
        return self._epsilon

    def decode(self, seq: bytes) -> DecodedCodon:
        ccode = ffi.new("struct nmm_codon *")
        lprob: float = lib.nmm_frame_state_decode(self._cdata, seq, len(seq), ccode)
        return DecodedCodon(lprob, ccode.a + ccode.b + ccode.c)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.nmm_frame_state_destroy(self._cdata)
