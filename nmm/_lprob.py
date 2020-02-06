from math import inf, nan, isnan, isinf

LPROB_ZERO: float = -inf
LPROB_INVALID: float = nan


def lprob_is_zero(x: float):
    return isinf(x) and x < 0


def lprob_is_valid(x: float):
    return not isnan(x)
