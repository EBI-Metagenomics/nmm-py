from ._frame_state import FrameState
from ._hmm import HMM
from ._log import LOG
from ._state import NormalState, MuteState, TableState
from ._alphabet import Alphabet

_ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

try:
    from ._ffi import ffi as _

    assert _ is not None
except Exception as e:
    e.msg = e.msg + _ffi_err
    raise e

__version__ = "0.0.1"

__all__ = [
    "__version__",
    "HMM",
    "MuteState",
    "NormalState",
    "TableState",
    "FrameState",
    "LOG",
    "Alphabet",
]
