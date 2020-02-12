from typing import Iterator, Tuple

from ..._subsequence import CSubSequence
from ..fragment import Fragment
from .path import StandardPath
from .step import StandardStep


class StandardFragment(Fragment):
    """
    Fragment of the standard profile.

    Parameters
    ----------
    subsequence : `CSubSequence`
        Sequence.
    path : `StandardPath`
        Path of the standard profile.
    homologous : `bool`
        Fragment homology.
    """

    def __init__(
        self, subsequence: CSubSequence, path: StandardPath, homologous: bool,
    ):
        super().__init__(homologous)
        self._subsequence = subsequence
        self._path = path

    @property
    def subsequence(self) -> CSubSequence:
        return self._subsequence

    def items(self) -> Iterator[Tuple[bytes, StandardStep]]:
        start = end = 0
        for step in self._path:
            end += step.seq_len
            yield (self._subsequence.symbols[start:end], step)
            start = end

    def __repr__(self):
        seq = self.sequence.decode()
        return f"<{self.__class__.__name__}:{seq}>"
