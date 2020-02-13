from math import inf, isinf, isnan, nan

LPROB_ZERO: float = -inf
LPROB_INVALID: float = nan


def lprob_is_zero(x: float):
    return isinf(x) and x < 0


def lprob_is_valid(x: float):
    return not isnan(x)


def lprob_normalize(arr):
    from numpy import asarray
    from scipy.special import logsumexp

    arr = asarray(arr, float)
    return arr - logsumexp(arr)
