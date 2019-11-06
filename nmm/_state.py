from typing import Dict

from ._alphabet import Alphabet
from ._base import Base
from ._codon import Codon
from ._ffi import ffi, lib
from ._log import LOG0


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
        seq : str
            Sequence.
        """
        return lib.imm_state_lprob(self.__cdata, seq, len(seq))


class State(CState):
    def __init__(self, cdata: ffi.CData, alphabet: Alphabet):
        """
        Parameters
        ----------
        alphabet : str
            Alphabet.
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
        name : bytes
            Name.
        alphabet : Alphabet
            Alphabet.
        """
        cdata = lib.imm_mute_state_create(name, alphabet.cdata)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")
        self._cdata = cdata
        super(MuteState, self).__init__(lib.imm_state_cast_c(cdata), alphabet)

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
        name : bytes
            Name.
        alphabet : Alphabet
            Alphabet.
        lprobs : dict
            List of probabilities in log-space.
        """

        if len(set(b"".join(lprobs.keys())) - set(alphabet.symbols)) > 0:
            raise ValueError("Unrecognized alphabet symbol.")

        arr = [lprobs.get(bytes([symb]), LOG0) for symb in alphabet.symbols]
        cdata = lib.imm_normal_state_create(name, alphabet.cdata, arr)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        self._cdata = cdata
        super(NormalState, self).__init__(lib.imm_state_cast_c(cdata), alphabet)

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
        name : bytes
            Name.
        alphabet : Alphabet
            Alphabet.
        emission : dict
            Emission probabilities in log-space.
        """
        cdata = lib.imm_table_state_create(name, alphabet.cdata)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_table_state_create` failed.")

        for seq, lprob in emission.items():
            lib.imm_table_state_add(cdata, seq, lprob)

        self._cdata = cdata
        super(TableState, self).__init__(lib.imm_state_cast_c(cdata), alphabet)

    def normalize(self) -> None:
        err = lib.imm_table_state_normalize(self._cdata)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.imm_table_state_destroy(self._cdata)


class FrameState(State):
    def __init__(self, name: bytes, base: Base, codon: Codon, epsilon: float):
        """
        Parameters
        ----------
        name : bytes
            Name.
        base : Base
            Base.
        codon : Codon
            Codon.
        epsilon : float
            Epsilon.
        """
        if set(base.alphabet.symbols) != set(codon.alphabet.symbols):
            raise ValueError("Alphabet symbols of `base` and `codon` are not equal.")

        self._base = base
        self._codon = codon
        self._epsilon = epsilon

        cdata = lib.nmm_frame_state_create(name, base.cdata, codon.cdata, epsilon)
        if cdata == ffi.NULL:
            raise RuntimeError("Could not create state.")

        self._cdata = cdata
        super(FrameState, self).__init__(lib.imm_state_cast_c(cdata), codon.alphabet)

    @property
    def base(self) -> Base:
        return self._base

    @property
    def codon(self) -> Codon:
        return self._codon

    @property
    def epsilon(self):
        return self._epsilon

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name.decode()}>"

    def __del__(self):
        if self._cdata != ffi.NULL:
            lib.nmm_frame_state_destroy(self._cdata)
