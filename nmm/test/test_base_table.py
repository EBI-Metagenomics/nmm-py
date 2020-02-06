from math import log, isnan
from nmm import Alphabet, Base, BaseTable

import pytest
from numpy.testing import assert_allclose, assert_equal


def test_base_table():
    base = Base(Alphabet(b"ACGT", b"X"))
    baset = BaseTable(base, (log(0.1), log(0.2), log(0.3), log(0.4)))
    assert_allclose(baset.lprob(b"A"), log(0.1))
    assert_allclose(baset.lprob(b"C"), log(0.2))
    assert_allclose(baset.lprob(b"G"), log(0.3))
    assert_allclose(baset.lprob(b"T"), log(0.4))
    assert_equal(isnan(baset.lprob(b"X")), True)

    with pytest.raises(Exception):
        baset = BaseTable(base, (log(0.1), log(0.2), log(0.3)))
