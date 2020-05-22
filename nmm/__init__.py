from ._alphabet import (
    AlphabetType,
    AminoAlphabet,
    BaseAlphabet,
    CanonicalAminoAlphabet,
    DNAAlphabet,
    RNAAlphabet,
)
from ._cdata import CData
from ._codon import Codon, codon_iter
from ._codon_prob import CodonProb
from ._gencode import GeneticCode
from ._input import Input
from ._model import Model
from ._output import Output
from ._state import CodonState, FrameState, StateType
from ._table import AminoTable, BaseTable, CodonTable
from ._testit import test

try:
    from ._ffi import lib
except Exception as e:
    _ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

    raise RuntimeError(str(e) + _ffi_err)

__version__ = "0.0.6"

__all__ = [
    "AlphabetType",
    "AminoAlphabet",
    "AminoTable",
    "BaseAlphabet",
    "BaseTable",
    "CData",
    "CanonicalAminoAlphabet",
    "Codon",
    "CodonProb",
    "CodonState",
    "CodonTable",
    "DNAAlphabet",
    "FrameState",
    "GeneticCode",
    "Input",
    "Model",
    "Output",
    "RNAAlphabet",
    "StateType",
    "__version__",
    "codon_iter",
    "lib",
    "test",
]
