from math import log

import pytest
from numpy import nan
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, BaseTable, Base


def test_base():
    alphabet = Alphabet(b"ACGT")
    base = BaseTable.create(alphabet)
    base.set_lprob(Base("A"), log(0.3))
    base.set_lprob(Base("T"), log(0.3))

    assert_allclose(base.get_lprob(Base("A")), log(0.3))
    assert_allclose(base.get_lprob(Base("T")), log(0.3))

    assert_equal(base.get_lprob(Base("C")), LOG0)
    assert_equal(base.get_lprob(Base("G")), LOG0)

    assert_equal(base.get_lprob(Base("X")), nan)

    with pytest.raises(ValueError):
        base.set_lprob(Base("X"), 0.0)

    base.normalize()

    assert_allclose(base.get_lprob(Base("A")), log(0.3) - log(0.6))
    assert_allclose(base.get_lprob(Base("T")), log(0.3) - log(0.6))

    assert_equal(base.get_lprob(Base("C")), LOG0)
    assert_equal(base.get_lprob(Base("G")), LOG0)

    assert_equal(base.get_lprob(Base("X")), nan)

    assert_equal(set(base.alphabet.symbols), set(b"ACGT"))

    base.set_lprob(Base("A"), LOG0)
    base.set_lprob(Base("C"), LOG0)
    base.set_lprob(Base("G"), LOG0)
    base.set_lprob(Base("T"), LOG0)

    with pytest.raises(RuntimeError):
        base.normalize()
