from typing import Iterator, List, Tuple, TypeVar, Union

from .._ffi import ffi
from .._gencode import GeneticCode
from .._log import LOG1
from .._path import Path
from .._state import CodonState, MuteState, State, TableState, NormalState
from .._step import Step
from .._alphabet import Alphabet
from .result import Fragment

from .amino_acid import AminoAcidFragment, AminoAcidPath


class CodonStep(Step):
    def __init__(
        self, imm_step: ffi.CData, state: Union[MuteState, CodonState], seq_len: int
    ):
        super().__init__(imm_step, state, seq_len)
        self._state = state

    @property
    def state(self) -> Union[MuteState, CodonState]:
        return self._state


T = TypeVar("T", bound="CodonPath")


class CodonPath(Path):
    def __init__(self):
        super().__init__()
        self._steps: List[CodonStep] = []

    def append_codon_step(
        self, state: Union[MuteState, CodonState], seq_len: int
    ) -> ffi.CData:
        imm_step = self._append_imm_step(state.imm_state, seq_len)
        step = CodonStep(imm_step, state, seq_len)
        self._steps.append(step)
        return step

    def append(self, state: State, seq_len: int) -> ffi.CData:
        # TODO: think in a better solution
        # Solution: I should not have a class base with methods
        # i cannot truly inherit.
        raise RuntimeError("Call `append_codon_step` instead.")
        del state
        del seq_len

    def steps(self) -> Iterator[CodonStep]:
        return iter(self._steps)


class CodonFragment(Fragment):
    def __init__(
        self, sequence: bytes, path: CodonPath, homologous: bool,
    ):
        super().__init__(sequence, homologous)
        self._path = path

    def items(self) -> Iterator[Tuple[bytes, CodonStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self.sequence[start:end], step)
            start = end

    def decode(self, genetic_code: GeneticCode) -> AminoAcidFragment:
        nseq: List[bytes] = []
        npath = AminoAcidPath()

        amino_acids = genetic_code.amino_acids()
        aa_alphabet = Alphabet(b"".join(amino_acids))

        start: int = 0
        seq = self.sequence
        for step in self._path.steps():
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, aa_alphabet)
                npath.append_amino_acid_step(mstate, 0)
            else:
                assert isinstance(step.state, TableState)
                aa = genetic_code.amino_acid(seq[start : start + step.seq_len])
                nseq.append(aa)

                nstate = NormalState(step.state.name, aa_alphabet, {aa: LOG1})
                npath.append_amino_acid_step(nstate, 1)

            start += step.seq_len

        return AminoAcidFragment(b"".join(nseq), npath, self.homologous)

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"
