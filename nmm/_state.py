from abc import ABC, abstractmethod
from typing import Dict, NamedTuple

from ._alphabet import Alphabet
from ._base import BaseTable
from ._codon import CodonTable
from ._ffi import ffi, lib
from ._log import LOG0

DecodedCodon = NamedTuple("DecodedCodon", [("lprob", float), ("codon", bytes)])


class State(ABC):
    def __init__(self, alphabet: Alphabet):
        """
        Parameters
        ----------
        cdata : `struct imm_state *`
        alphabet : Alphabet.
        """
        self._alphabet = alphabet

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    @property
    @abstractmethod
    def name(self) -> bytes:
        raise NotImplementedError()

    @property
    @abstractmethod
    def min_seq(self) -> int:
        raise NotImplementedError()

    @property
    @abstractmethod
    def max_seq(self) -> int:
        raise NotImplementedError()

    @abstractmethod
    def lprob(self, seq: bytes) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        seq : Sequence.
        """
        del seq
        raise NotImplementedError()

    def __str__(self) -> str:
        return f"<{self.name.decode()}>"


class CState:
    def __init__(self, imm_state: ffi.CData):
        self.__imm_state = imm_state

    @property
    def name(self) -> bytes:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self.__imm_state == ffi.NULL:
            raise RuntimeError("State has failed to initialize.")
        return ffi.string(lib.imm_state_get_name(self.__imm_state))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self.__imm_state)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self.__imm_state)

    @property
    def imm_state(self) -> ffi.CData:
        return self.__imm_state

    def lprob(self, seq: bytes) -> float:
        return lib.imm_state_lprob(self.__imm_state, seq, len(seq))


class MuteState(CState, State):
    def __init__(self, name: bytes, alphabet: Alphabet):
        """
        Parameters
        ----------
        name : Name.
        alphabet : Alphabet.
        """
        self._imm_state = ffi.NULL
        self._imm_state = lib.imm_mute_state_create(name, alphabet.imm_abc)
        if self._imm_state == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")

        CState.__init__(self, lib.imm_state_cast_c(self._imm_state))
        State.__init__(self, alphabet)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._imm_state != ffi.NULL:
            lib.imm_mute_state_destroy(self._imm_state)


class NormalState(CState, State):
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

        self._imm_state = ffi.NULL
        self._imm_state = lib.imm_normal_state_create(name, alphabet.imm_abc, arr)
        if self._imm_state == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        CState.__init__(self, lib.imm_state_cast_c(self._imm_state))
        State.__init__(self, alphabet)

    def emission_table(self) -> Dict[bytes, float]:
        return {bytes([s]): self.lprob(bytes([s])) for s in self.alphabet.symbols}

    def normalize(self) -> None:
        err = lib.imm_normal_state_normalize(self._imm_state)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._imm_state != ffi.NULL:
            lib.imm_normal_state_destroy(self._imm_state)


class TableState(CState, State):
    def __init__(self, name: bytes, alphabet: Alphabet, emission: Dict[bytes, float]):
        """
        Parameters
        ----------
        name : Name.
        alphabet : Alphabet.
        emission : Emission probabilities in log-space.
        """
        self._imm_state = ffi.NULL
        self._imm_state = lib.imm_table_state_create(name, alphabet.imm_abc)
        if self._imm_state == ffi.NULL:
            raise RuntimeError("`imm_table_state_create` failed.")

        for seq, lprob in emission.items():
            lib.imm_table_state_add(self._imm_state, seq, lprob)

        CState.__init__(self, lib.imm_state_cast_c(self._imm_state))
        State.__init__(self, alphabet)

    def normalize(self) -> None:
        err = lib.imm_table_state_normalize(self._imm_state)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._imm_state != ffi.NULL:
            lib.imm_table_state_destroy(self._imm_state)


class CodonState(TableState):
    # TODO: consider creating a Codon type and place it in the keys as a type.
    def __init__(self, name: bytes, alphabet: Alphabet, emission: Dict[bytes, float]):
        if sum(len(k) != 3 for k in emission.keys()) > 0:
            raise ValueError("Codon must be composed of three bases.")

        super().__init__(name, alphabet, emission)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"


class FrameState(CState, State):
    def __init__(
        self, name: bytes, baset: BaseTable, codont: CodonTable, epsilon: float
    ):
        """
        Parameters
        ----------
        name : Name.
        baset : Base table.
        codont : Codon table.
        epsilon : Epsilon.
        """
        if set(baset.alphabet.symbols) != set(codont.alphabet.symbols):
            raise ValueError("Alphabet symbols of `base` and `codon` are not equal.")

        self._baset = baset
        self._codont = codont
        self._epsilon = epsilon

        self._imm_state = ffi.NULL
        self._imm_state = lib.nmm_frame_state_create(
            name, baset.nmm_baset, codont.nmm_codont, epsilon
        )
        if self._imm_state == ffi.NULL:
            raise RuntimeError("Could not create state.")

        CState.__init__(self, lib.imm_state_cast_c(self._imm_state))
        State.__init__(self, codont.alphabet)

    @property
    def base(self) -> BaseTable:
        return self._baset

    @property
    def codon(self) -> CodonTable:
        return self._codont

    @property
    def epsilon(self):
        return self._epsilon

    def decode(self, seq: bytes) -> DecodedCodon:
        ccode = ffi.new("struct nmm_codon *")
        lprob: float = lib.nmm_frame_state_decode(self._imm_state, seq, len(seq), ccode)
        return DecodedCodon(lprob, ccode.a + ccode.b + ccode.c)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._imm_state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._imm_state)
