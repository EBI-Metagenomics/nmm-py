from typing import Iterator, List, Tuple, TypeVar, Union

from .result import Fragment
from .._state import NormalState, MuteState, State

from .._ffi import ffi
from .._step import Step
from .._path import Path


class AminoAcidStep(Step):
    def __init__(
        self, imm_step: ffi.CData, state: Union[MuteState, NormalState], seq_len: int
    ):
        super().__init__(imm_step, state, seq_len)
        self._state = state

    @property
    def state(self) -> Union[MuteState, NormalState]:
        return self._state


T = TypeVar("T", bound="AminoAcidPath")


class AminoAcidPath(Path):
    def __init__(self):
        super().__init__()
        self._steps: List[AminoAcidStep] = []

    def append_amino_acid_step(
        self, state: Union[MuteState, NormalState], seq_len: int
    ) -> ffi.CData:
        imm_step = self._append_imm_step(state.imm_state, seq_len)
        step = AminoAcidStep(imm_step, state, seq_len)
        self._steps.append(step)
        return step

    def append(self, state: State, seq_len: int) -> ffi.CData:
        # TODO: think in a better solution
        # Solution: I should not have a class base with methods
        # i cannot truly inherit.
        raise RuntimeError("Call `append_amino_acid_step` instead.")
        del state
        del seq_len

    def steps(self) -> Iterator[AminoAcidStep]:
        return iter(self._steps)


class AminoAcidFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: AminoAcidPath, homologous: bool,
    ):
        super().__init__(sequence, homologous)
        self._path = path

    def items(self) -> Iterator[Tuple[bytes, AminoAcidStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"
