from . import alphabet, codon, fragment, path, prob, result, sequence, state
from ._gencode import GeneticCode
from ._imm import (
    HMM,
    wrap_imm_path,
    wrap_imm_result,
    wrap_imm_results,
)
from ._interval import Interval
from ._testit import test

try:
    from ._ffi import ffi

    CData = ffi.CData
    del ffi
except Exception as e:
    _ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

    raise RuntimeError(str(e) + _ffi_err)

__version__ = "0.0.4"

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
    "wrap_imm_path",
    "wrap_imm_result",
    "wrap_imm_results",
]
