import pytest
from numpy import nan
from numpy.testing import assert_allclose, assert_equal
from nmm import Base, Alphabet, LOG


def test_base():
    alphabet = Alphabet("ACGT")
    base = Base(alphabet)
    base.set_lprob("A", LOG(0.3))
    base.set_lprob("T", LOG(0.3))

    assert_allclose(base.get_lprob("A"), LOG(0.3))
    assert_allclose(base.get_lprob("T"), LOG(0.3))

    assert_equal(base.get_lprob("C"), LOG(0.0))
    assert_equal(base.get_lprob("G"), LOG(0.0))

    assert_equal(base.get_lprob("X"), nan)

    with pytest.raises(ValueError):
        base.set_lprob("X", 0.0)

    with pytest.raises(ValueError):
        base.set_lprob("XX", 0.0)

    with pytest.raises(ValueError):
        base.get_lprob("XX")

    base.normalize()

    assert_allclose(base.get_lprob("A"), LOG(0.3) - LOG(0.6))
    assert_allclose(base.get_lprob("T"), LOG(0.3) - LOG(0.6))

    assert_equal(base.get_lprob("C"), LOG(0.0))
    assert_equal(base.get_lprob("G"), LOG(0.0))

    assert_equal(base.get_lprob("X"), nan)
