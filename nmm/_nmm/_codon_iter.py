from typing import Iterable
import itertools

from ._base_alphabet import BaseAlphabet
from ._codon import Codon


def codon_iter(base_abc: BaseAlphabet) -> Iterable[Codon]:
    """
    Codon iterator.

    Parameters
    ----------
    base_abc
        Base alphabet.
    """
    bases = [base_abc.symbols[i : i + 1] for i in range(len(base_abc.symbols))]

    for a, b, c in itertools.product(bases, bases, bases):
        yield Codon(a + b + c, base_abc)
