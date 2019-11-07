from ._cli import cli
from ._alphabet import Alphabet
from ._base import Base
from ._codon import Codon
from ._hmm import HMM
from ._hmmer import create_frame_profile, create_hmmer_profile, read_hmmer
from ._log import LOG0
from ._path import Path
from ._state import FrameState, MuteState, NormalState, TableState
from ._step import Step
from ._fasta import FASTAReader, FASTAParserError

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
    "FASTAReader",
    "FASTAParserError",
    "Alphabet",
    "Base",
    "Codon",
    "FrameState",
    "HMM",
    "LOG0",
    "MuteState",
    "NormalState",
    "Path",
    "Step",
    "TableState",
    "__version__",
    "create_frame_profile",
    "create_hmmer_profile",
    "read_hmmer",
    "cli",
]
