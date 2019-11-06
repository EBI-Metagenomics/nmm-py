from typing import List
from .._step import CStep


class Fragment:
    def __init__(self, seq: bytes, steps: List[CStep] = []):
        self._seq = seq
        self._steps = steps

    @property
    def sequence(self) -> bytes:
        return self._seq


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
