from typing import List, Tuple, Iterator
from .._step import CStep
from .._path import CPath


class Fragment:
    def __init__(self, seq: bytes, steps: List[CStep] = []):
        self._seq = seq
        self._path = CPath(None)
        for step in steps:
            self._path.append(step.state.imm_state, step.seq_len)

    @property
    def sequence(self) -> bytes:
        return self._seq

    def items(self) -> Iterator[Tuple[bytes, CStep]]:
        start = end = 0
        for step in self._path.steps():
            end += step.seq_len
            yield (self._seq[start:end], step)
            start = end


class HomoFragment(Fragment):
    def __init__(self, seq: bytes, steps: List[CStep] = []):
        super().__init__(seq, steps)

    @property
    def homologous(self):
        return True

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"


class NonHomoFragment(Fragment):
    def __init__(self, seq: bytes, steps: List[CStep] = []):
        super().__init__(seq, steps)

    @property
    def homologous(self):
        return False

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"
