from ._alphabet import Alphabet
from ._log import LOG
from ._string import make_sure_bytes

from ._ffi import ffi, lib


class State:
    def __init__(self, alphabet: Alphabet):
        """
        Parameters
        ----------
        name : str
            Name.
        alphabet : str
            Alphabet.
        """
        self._state = ffi.NULL
        self._alphabet = alphabet

    @property
    def name(self) -> str:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self._state == ffi.NULL:
            return "unknown-name"
        n = ffi.string(lib.imm_state_get_name(self.cdata))
        return n.decode()

    @property
    def alphabet(self) -> Alphabet:
        return self._alphabet

    def lprob(self, seq: str) -> float:
        """
        Log-probability of sequence emission.

        Parameters
        ----------
        seq : str
            Sequence.
        """
        return lib.imm_state_lprob(self.cdata, make_sure_bytes(seq), len(seq))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self.cdata)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self.cdata)

    @property
    def cdata(self):
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
        name = make_sure_bytes(name)
        self._state = lib.imm_mute_state_create(name, alphabet.cdata)

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_mute_state_destroy(self._state)

    def normalize(self):
        pass

    def emission(self):
        return [("", LOG(1.0))]

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"


class NormalState(State):
    def __init__(self, name: str, alphabet: Alphabet, lprobs: dict):
        """
        Parameters
        ----------
        name : str
            Name.
        alphabet : Alphabet
            Alphabet.
        lprobs : list
            List of probabilities in log-space.
        """
        super(NormalState, self).__init__(alphabet)

        if len(set(lprobs.keys() - set(alphabet.symbols))) > 0:
            raise ValueError("Unrecognized alphabet symbol.")

        lprobs_arr = [lprobs.get(symb, LOG(0.0)) for symb in alphabet.symbols]

        name = make_sure_bytes(name)
        self._state = lib.imm_normal_state_create(name, alphabet.cdata, lprobs_arr)

    def __del__(self):
        if self._state != ffi.NULL:
            lib.imm_normal_state_destroy(self._state)

    def normalize(self):
        err = lib.imm_normal_state_normalize(self._state)
        if err != 0:
            raise ValueError("Normalization error.")

    def emission(self):
        emission = {symbol: self.lprob(symbol) for symbol in self.alphabet.symbols}
        return emission_table(emission)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"


class TableState(State):
    def __init__(self, name: str, alphabet: Alphabet, emission: dict):
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

        name = make_sure_bytes(name)
        self._state = lib.imm_table_state_create(name, alphabet.cdata)

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


def emission_table(emission: dict):
    """
    Parameters
    ----------
    emission : dict
        Emission probabilities in log-space.
    """
    table = list(emission.items())
    table = sorted(table, key=lambda x: -x[1])
    return table
