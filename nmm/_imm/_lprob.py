from math import inf, isinf, isnan, nan
from typing import Iterable


def lprob_zero() -> float:
    return -inf


def lprob_invalid() -> float:
    return nan


def lprob_is_zero(x: float):
    return isinf(x) and x < 0


def lprob_is_valid(x: float):
    return not isnan(x)


def lprob_normalize(arr: Iterable[float]):
    from numpy import asarray
    from scipy.special import logsumexp

    arr = asarray(arr, float)
    return arr - logsumexp(arr)
