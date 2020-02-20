from . import alphabet, codon, fragment, path, prob, result, sequence, state
from ._gencode import GeneticCode
from ._imm import HMM
from ._interval import Interval
from ._testit import test
from ._cdata import CData

try:
    from ._ffi import ffi

    del ffi
except Exception as e:
    _ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

    raise RuntimeError(str(e) + _ffi_err)

__version__ = "0.0.5"

__all__ = [
    "GeneticCode",
    "HMM",
    "Interval",
    "__version__",
    "alphabet",
    "codon",
    "fragment",
    "path",
    "prob",
    "result",
    "sequence",
    "state",
    "test",
    "CData",
]
