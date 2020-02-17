from typing import Iterable
import itertools

from ._base_alphabet import CBaseAlphabet
from ._codon import Codon


def codon_iter(base_abc: CBaseAlphabet) -> Iterable[Codon]:
    """
    Codon iterator.

    Parameters
    ----------
    base_abc : CBaseAlphabet
        Base alphabet.
    """
    bases = [base_abc.symbols[i : i + 1] for i in range(len(base_abc.symbols))]

    for a, b, c in itertools.product(bases, bases, bases):
        yield Codon(a + b + c, base_abc)
