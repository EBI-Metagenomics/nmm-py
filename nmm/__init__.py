from ._alphabet import Alphabet
from ._base import Base
from ._codon import Codon
from ._sequence import Sequence
from ._sequence_table import SequenceTable
from ._base_table import BaseTable
from ._codon_prob import CodonProb
from ._codon_table import CodonTable
from ._lprob import lprob_is_zero, lprob_is_valid, LPROB_ZERO, LPROB_INVALID
from ._state import MuteState, NormalState, TableState

# from ._state import FrameState, MuteState, NormalState, TableState, CodonState

# from ._base import BaseTable, Base
# from ._cli import cli
# from ._codon import CodonTable, Codon
# from ._gencode import GeneticCode
# from ._hmm import HMM
# from ._hmmer import create_frame_profile, create_standard_profile
# from ._hmmer.io import tblout_reader
# from ._log import LOG0
# from ._path import CPath
# from ._step import CStep
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

__all__ = [
    "Alphabet",
    "Base",
    "BaseTable",
    "Codon",
    "CodonProb",
    "CodonTable",
    "LPROB_INVALID",
    "LPROB_ZERO",
    "Sequence",
    "__version__",
    "cli",
    "lprob_is_valid",
    "lprob_is_zero",
    "test",
    "FrameState",
    "MuteState",
    "NormalState",
    "TableState",
    "CodonState",
    "SequenceTable",
    "TableState",
]
