from typing import Sequence

from ._sequence import CSequence
from ._alphabet import CAlphabet
from ._lprob import lprob_is_valid

# from ._base import BaseTable
# from ._codon import Codon, CodonTable
from ._ffi import ffi, lib

# from ._log import LOG0


class CState:
    def __init__(self, imm_state: ffi.CData):
        if imm_state == ffi.NULL:
            raise RuntimeError("`imm_state` is NULL.")
        self._imm_state = imm_state

    @property
    def imm_state(self) -> ffi.CData:
        return self._imm_state

    @property
    def name(self) -> bytes:
        return ffi.string(lib.imm_state_get_name(self._imm_state))

    @property
    def min_seq(self) -> int:
        return lib.imm_state_min_seq(self._imm_state)

    @property
    def max_seq(self) -> int:
        return lib.imm_state_max_seq(self._imm_state)

    def lprob(self, seq: CSequence) -> float:
        """
        Log-space probability of sequence emission.

        Parameters
        ----------
        seq : `CSequence`
            Sequence.
        """
        lprob: float = lib.imm_state_lprob(self._imm_state, seq.imm_seq)
        if not lprob_is_valid(lprob):
            raise RuntimeError("Could not get probability.")
        return lprob

    def __str__(self) -> str:
        # Refer to https://github.com/pytest-dev/pytest/issues/4659
        if self._imm_state == ffi.NULL:
            raise RuntimeError("State has failed to initialize.")
        return f"{self.name.decode()}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class MuteState(CState):
    def __init__(self, name: bytes, alphabet: CAlphabet):
        """
        Mute state.

        Parameters
        ----------
        name : bytes
            State name.
        alphabet : `CAlphabet`
            Alphabet.
        """
        self._imm_mute_state = lib.imm_mute_state_create(name, alphabet.imm_abc)
        if self._imm_mute_state == ffi.NULL:
            raise RuntimeError("`imm_mute_state_create` failed.")

        super().__init__(lib.imm_state_cast_c(self._imm_mute_state))
        self._alphabet = alphabet

    @property
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    def __del__(self):
        if self._imm_mute_state != ffi.NULL:
            lib.imm_mute_state_destroy(self._imm_mute_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


class NormalState(CState):
    def __init__(self, name: bytes, alphabet: CAlphabet, lprobs: Sequence[float]):
        """
        Normal state.

        Parameters
        ----------
        name : bytes
            State name.
        alphabet : `CAlphabet`
            Alphabet.
        lprobs : `typing.Sequence[float]`
            Emission probabilities in log-space for each alphabet letter.
        """
        state = lib.imm_normal_state_create(name, alphabet.imm_abc, list(lprobs))
        if state == ffi.NULL:
            raise RuntimeError("`imm_normal_state_create` failed.")

        self._imm_normal_state = state
        super().__init__(lib.imm_state_cast_c(self._imm_normal_state))
        self._alphabet = alphabet

    @property
    def alphabet(self) -> CAlphabet:
        return self._alphabet

    def __del__(self):
        if self._imm_normal_state != ffi.NULL:
            lib.imm_normal_state_destroy(self._imm_normal_state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{str(self)}>"


# class TableState(CState):
#     def __init__(self, name: bytes, alphabet: CAlphabet, emission: Dict[bytes, float]):
#         """
#         Parameters
#         ----------
#         name : Name.
#         alphabet : CAlphabet.
#         emission : Emission probabilities in log-space.
#         """
#         self._imm_table_state = lib.imm_table_state_create(name, alphabet.imm_abc)
#         if self._imm_table_state == ffi.NULL:
#             raise RuntimeError("`imm_table_state_create` failed.")

#         for seq, lprob in emission.items():
#             lib.imm_table_state_add(self._imm_table_state, seq, lprob)

#         super().__init__(lib.imm_state_cast_c(self._imm_table_state))
#         self._alphabet = alphabet

#     @property
#     def alphabet(self) -> CAlphabet:
#         return self._alphabet

#     def normalize(self) -> None:
#         err = lib.imm_table_state_normalize(self._imm_table_state)
#         if err != 0:
#             raise RuntimeError("Normalization error.")

#     def __repr__(self):
#         return f"<{self.__class__.__name__}:{self.name.decode()}>"

#     def __del__(self):
#         if self._imm_table_state != ffi.NULL:
#             lib.imm_table_state_destroy(self._imm_table_state)


# class CodonState(TableState):
#     def __init__(self, name: bytes, alphabet: CAlphabet, emission: Dict[Codon, float]):
#         super().__init__(name, alphabet, {bytes(k): v for k, v in emission.items()})

#     def __repr__(self):
#         return f"<{self.__class__.__name__}:{self.name.decode()}>"


# class FrameState(CState):
#     def __init__(
#         self, name: bytes, baset: BaseTable, codont: CodonTable, epsilon: float
#     ):
#         """
#         Parameters
#         ----------
#         name : Name.
#         baset : Base table.
#         codont : Codon table.
#         epsilon : Epsilon.
#         """
#         if set(baset.alphabet.symbols) != set(codont.alphabet.symbols):
#             raise ValueError("Alphabet symbols of `base` and `codon` are not equal.")

#         self._baset = baset
#         self._codont = codont
#         self._epsilon = epsilon

#         self._imm_frame_state = lib.nmm_frame_state_create(
#             name, baset.nmm_baset, codont.nmm_codont, epsilon
#         )
#         if self._imm_frame_state == ffi.NULL:
#             raise RuntimeError("Could not create state.")

#         super().__init__(lib.imm_state_cast_c(self._imm_frame_state))

#     @property
#     def alphabet(self) -> CAlphabet:
#         return self._baset.alphabet

#     @property
#     def base(self) -> BaseTable:
#         return self._baset

#     @property
#     def codon(self) -> CodonTable:
#         return self._codont

#     @property
#     def epsilon(self):
#         return self._epsilon

#     def decode(self, seq: bytes) -> Tuple[Codon, float]:
#         codon = Codon("XXX")
#         st = self._imm_frame_state
#         lprob: float = lib.nmm_frame_state_decode(st, seq, len(seq), codon.nmm_codon)
#         return (codon, lprob)

#     def __repr__(self):
#         return f"<{self.__class__.__name__}:{self.name.decode()}>"

#     def __del__(self):
#         if self._imm_frame_state != ffi.NULL:
#             lib.nmm_frame_state_destroy(self._imm_frame_state)
