from math import log
import pytest
from numpy import nan
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, Base


def test_base():
    alphabet = Alphabet("ACGT")
    base = Base(alphabet)
    base.set_lprob("A", log(0.3))
    base.set_lprob("T", log(0.3))

    assert_allclose(base.get_lprob("A"), log(0.3))
    assert_allclose(base.get_lprob("T"), log(0.3))

    assert_equal(base.get_lprob("C"), LOG0)
    assert_equal(base.get_lprob("G"), LOG0)

    assert_equal(base.get_lprob("X"), nan)

    with pytest.raises(ValueError):
        base.set_lprob("X", 0.0)

    with pytest.raises(ValueError):
        base.set_lprob("XX", 0.0)

    with pytest.raises(ValueError):
        base.get_lprob("XX")

    base.normalize()

    assert_allclose(base.get_lprob("A"), log(0.3) - log(0.6))
    assert_allclose(base.get_lprob("T"), log(0.3) - log(0.6))

    assert_equal(base.get_lprob("C"), LOG0)
    assert_equal(base.get_lprob("G"), LOG0)

    assert_equal(base.get_lprob("X"), nan)

    assert_equal(set(base.alphabet.symbols), set("ACGT"))

    base.set_lprob("A", LOG0)
    base.set_lprob("C", LOG0)
    base.set_lprob("G", LOG0)
    base.set_lprob("T", LOG0)

    with pytest.raises(RuntimeError):
        base.normalize()
