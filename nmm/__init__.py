from ._alphabet import Alphabet
from ._base import Base
from ._codon import Codon
from ._hmm import HMM
from ._hmmer import create_frame_profile, create_hmmer_profile, read_hmmer
from ._log import LOG0
from ._path import Path
from ._state import FrameState, MuteState, NormalState, TableState

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
    "__version__",
    "HMM",
    "MuteState",
    "NormalState",
    "TableState",
    "FrameState",
    "Alphabet",
    "Codon",
    "Base",
    "Path",
    "LOG0",
    "read_hmmer",
    "create_frame_profile",
    "create_hmmer_profile",
]
