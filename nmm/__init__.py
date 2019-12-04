from ._alphabet import Alphabet
from ._base import BaseTable, Base
from ._cli import cli
from ._codon import CodonTable, Codon
from ._gencode import GeneticCode
from ._hmm import HMM
from ._hmmer import create_frame_profile, create_standard_profile, read_hmmer
from ._log import LOG0
from ._path import CPath
from ._state import FrameState, MuteState, NormalState, TableState, CodonState
from ._step import CStep
from ._testit import test

try:
    from ._ffi import ffi as _

    del _
except Exception as e:
    _ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

    raise RuntimeError(str(e) + _ffi_err)

__version__ = "0.0.1"

__all__ = [
    "Alphabet",
    "Base",
    "BaseTable",
    "CPath",
    "CStep",
    "Codon",
    "CodonState",
    "CodonTable",
    "FrameState",
    "GeneticCode",
    "HMM",
    "LOG0",
    "MuteState",
    "NormalState",
    "TableState",
    "__version__",
    "cli",
    "create_frame_profile",
    "create_standard_profile",
    "read_hmmer",
    "test",
]
