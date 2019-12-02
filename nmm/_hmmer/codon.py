from typing import Iterator, List, Sequence, Tuple, TypeVar, Union

from .._alphabet import Alphabet
from .._codon import Codon
from .._ffi import ffi
from .._gencode import GeneticCode
from .._log import LOG1
from .._path import CPath
from .._state import CodonState, MuteState, NormalState, TableState
from .._step import CStep
from .amino_acid import AminoAcidFragment, AminoAcidPath
from .result import Fragment


class CodonStep(CStep):
    def __init__(self, imm_step: ffi.CData, state: Union[MuteState, CodonState]):
        super().__init__(imm_step)
        self._state = state

    @property
    def state(self) -> Union[MuteState, CodonState]:
        return self._state


T = TypeVar("T", bound="CodonPath")


class CodonPath(CPath):
    def __init__(self):
        super().__init__()
        self._steps: List[CodonStep] = []

    def append_codon_step(
        self, state: Union[MuteState, CodonState], seq_len: int
    ) -> CodonStep:
        cstep = self.append_cstep(state, seq_len)
        self._steps.append(CodonStep(cstep.imm_step, state))
        return self._steps[-1]

    def steps(self) -> Iterator[CodonStep]:
        return iter(self._steps)


class CodonFragment(Fragment):
    def __init__(
        self, codons: Sequence[Codon], path: CodonPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._codons = codons
        self._path = path

    @property
    def sequence(self) -> bytes:
        return b"".join(bytes(codon) for codon in self._codons)

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
        i = 0
        for step in self._path.steps():
            if isinstance(step.state, MuteState):
                mstate = MuteState(step.state.name, aa_alphabet)
                npath.append_amino_acid_step(mstate, 0)
            else:
                assert isinstance(step.state, TableState)
                aa = genetic_code.amino_acid(self._codons[i])
                nseq.append(aa)

                nstate = NormalState(step.state.name, aa_alphabet, {aa: LOG1})
                npath.append_amino_acid_step(nstate, 1)
                i += 1

            start += step.seq_len

        return AminoAcidFragment(b"".join(nseq), npath, self.homologous)

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"
