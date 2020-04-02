from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm.alphabet import BaseAlphabet
from nmm.prob import BaseTable, lprob_is_valid


def test_base_table():
    base = BaseAlphabet.create(b"ACGT", b"X")
    baset = BaseTable.create(base, (log(0.1), log(0.2), log(0.3), log(0.4)))
    assert_allclose(baset.lprob(b"A"), log(0.1))
    assert_allclose(baset.lprob(b"C"), log(0.2))
    assert_allclose(baset.lprob(b"G"), log(0.3))
    assert_allclose(baset.lprob(b"T"), log(0.4))

    with pytest.raises(Exception):
        baset = BaseTable.create(base, (log(0.1), log(0.2), log(0.3)))
