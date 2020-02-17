from ._gencode import GeneticCode
from ._imm import (
    HMM,
    LPROB_INVALID,
    LPROB_ZERO,
    Alphabet,
    AlphabetTable,
    CAlphabet,
    CAlphabetTable,
    CPath,
    CResult,
    CResults,
    CSequence,
    CSequenceTable,
    CState,
    CStep,
    CSubSequence,
    MuteState,
    NormalState,
    Path,
    Sequence,
    SequenceABC,
    SequenceTable,
    Step,
    SubSequence,
    TableState,
    create_imm_path,
    create_imm_step,
    lprob_is_valid,
    lprob_is_zero,
    lprob_normalize,
    wrap_imm_path,
    wrap_imm_result,
    wrap_imm_results,
)
from ._interval import Interval
from ._nmm import (
    AminoAlphabet,
    AminoTable,
    BaseAlphabet,
    BaseTable,
    CAminoAlphabet,
    CBaseAlphabet,
    Codon,
    CodonProb,
    CodonState,
    CodonTable,
    FrameState,
    codon_iter,
)
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
    "Alphabet",
    "AlphabetTable",
    "AminoAlphabet",
    "AminoTable",
    "BaseAlphabet",
    "BaseTable",
    "CAlphabet",
    "CAlphabetTable",
    "CAminoAlphabet",
    "CBaseAlphabet",
    "CData",
    "CPath",
    "CResult",
    "CResults",
    "CSequence",
    "CSequenceTable",
    "CState",
    "CStep",
    "CSubSequence",
    "Codon",
    "CodonProb",
    "CodonState",
    "CodonTable",
    "FrameState",
    "GeneticCode",
    "HMM",
    "Interval",
    "LPROB_INVALID",
    "LPROB_ZERO",
    "MuteState",
    "NormalState",
    "Path",
    "Sequence",
    "SequenceABC",
    "SequenceTable",
    "Step",
    "SubSequence",
    "TableState",
    "__version__",
    "codon_iter",
    "create_imm_path",
    "create_imm_step",
    "lprob_is_valid",
    "lprob_is_zero",
    "lprob_normalize",
    "test",
    "wrap_imm_path",
    "wrap_imm_result",
    "wrap_imm_results",
]
