from ._alphabet import Alphabet
from ._base import Base
from ._base_table import BaseTable
from ._codon import Codon
from ._codon_prob import CodonProb
from ._codon_table import CodonTable
from ._gencode import GeneticCode
from ._hmm import HMM

# from ._hmmer import create_frame_profile, create_standard_profile
# from ._hmmer.io import tblout_reader
from ._lprob import (
    LPROB_INVALID,
    LPROB_ZERO,
    lprob_is_valid,
    lprob_is_zero,
)

from ._path import Path
from ._results import CResult, CResults
from ._sequence import Sequence
from ._sequence_table import SequenceTable
from ._state import FrameState, MuteState, NormalState, TableState
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

__version__ = "0.0.4"

# from ._cli import cli

__all__ = [
    "Alphabet",
    "Base",
    "BaseTable",
    "CResult",
    "CResults",
    "Codon",
    "CodonProb",
    "CodonState",
    "CodonTable",
    "FrameState",
    "HMM",
    "LPROB_INVALID",
    "LPROB_ZERO",
    "MuteState",
    "NormalState",
    "Path",
    "Sequence",
    "SequenceTable",
    "TableState",
    "TableState",
    "__version__",
    "cli",
    # "create_frame_profile",
    # "create_standard_profile",
    "lprob_is_valid",
    "lprob_is_zero",
    "test",
    # "tblout_reader",
    "GeneticCode",
    "lprob_normalize",
]
