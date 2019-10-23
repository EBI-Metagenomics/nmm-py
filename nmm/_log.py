from math import inf, log

LOG0: float = -inf


def LOG(probability: float):
    if probability == 1.0:
        return 0.0
    if probability == 0.0:
        return -inf

    return log(probability)
