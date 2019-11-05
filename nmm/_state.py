from ._alphabet import Alphabet
from ._base import Base
from ._codon import Codon
from ._log import LOG0
from typing import Dict

from ._ffi import ffi, lib


class State:
    def __init__(self, alphabet: Alphabet):
        """
        Parameters
        ----------
        alphabet : str
            Alphabet.
        """
        self._state = ffi.NULL
        self._alphabet = alphabet

    @property
    def name(self) -> str:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self._state == ffi.NULL:
            raise RuntimeError("State has failed to initialize.")
        return ffi.string(lib.imm_state_get_name(self.cdata)).decode()

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    def lprob(self, seq: str) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        seq : str
            Sequence.
        """
        return lib.imm_state_lprob(self.cdata, seq.encode(), len(seq))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self.cdata)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self.cdata)

    @property
    def cdata(self) -> ffi.CData:
        return lib.imm_state_cast_c(self._state)

    def __str__(self) -> str:
        return f"<{self.name}>"


class MuteState(State):
    def __init__(self, name: str, alphabet: Alphabet):
        """
        Parameters
        ----------
        name : str
            Name.
        alphabet : Alphabet
            Alphabet.
        """
        super(MuteState, self).__init__(alphabet)
        self._state = lib.imm_mute_state_create(name.encode(), alphabet.cdata)
        if self._state == ffi.NULL:
            raise RuntimeError("Could not create state.")

    def normalize(self) -> None:
        pass

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_mute_state_destroy(self._state)


class NormalState(State):
    def __init__(self, name: str, alphabet: Alphabet, lprobs: Dict[str, float]):
        """
        Parameters
        ----------
        name : str
            Name.
        alphabet : Alphabet
            Alphabet.
        lprobs : dict
            List of probabilities in log-space.
        """
        super(NormalState, self).__init__(alphabet)

        if len(set(lprobs.keys() - set(alphabet.symbols))) > 0:
            raise ValueError("Unrecognized alphabet symbol.")

        arr = [lprobs.get(symb, LOG0) for symb in alphabet.symbols]
        self._state = lib.imm_normal_state_create(name.encode(), alphabet.cdata, arr)
        if self._state == ffi.NULL:
            raise RuntimeError("Could not create state.")

    def emission_table(self) -> Dict[str, float]:
        return {s: self.lprob(s) for s in self.alphabet.symbols}

    def normalize(self) -> None:
        err = lib.imm_normal_state_normalize(self._state)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_normal_state_destroy(self._state)


class TableState(State):
    def __init__(self, name: str, alphabet: Alphabet, emission: Dict[str, float]):
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
        super(TableState, self).__init__(alphabet)

        self._state = lib.imm_table_state_create(name.encode(), alphabet.cdata)
        if self._state == ffi.NULL:
            raise RuntimeError("Could not create state.")

        for seq, lprob in emission.items():
            lib.imm_table_state_add(self._state, seq.encode(), lprob)

    def normalize(self) -> None:
        err = lib.imm_table_state_normalize(self._state)
        if err != 0:
            raise RuntimeError("Normalization error.")

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_table_state_destroy(self._state)


class FrameState(State):
    def __init__(self, name: str, base: Base, codon: Codon, epsilon: float):
        """
        Parameters
        ----------
        name : str
            Name.
        base : Base
            Base.
        codon : Codon
            Codon.
        epsilon : float
            Epsilon.
        """
        super(FrameState, self).__init__(codon.alphabet)
        if set(base.alphabet.symbols) != set(codon.alphabet.symbols):
            raise ValueError("Alphabet symbols of `base` and `codon` are not equal.")

        self._base = base
        self._codon = codon
        self._epsilon = epsilon

        n = name.encode()
        self._state = lib.nmm_frame_state_create(n, base.cdata, codon.cdata, epsilon)
        if self._state == ffi.NULL:
            raise RuntimeError("Could not create state.")

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
        return f"<{self.__class__.__name__}:{self.name}>"

    def __del__(self):
        if self._state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._state)
