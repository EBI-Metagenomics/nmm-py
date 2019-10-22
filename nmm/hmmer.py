from pathlib import Path

import hmmer_reader

from ._log import LOG


def read_file(file):

    if not hasattr(file, "read"):
        if not isinstance(file, Path):
            file = Path(file)

        if not file.exists():
            raise ValueError(f"`{file}` does not exist.")

        if not file.is_file():
            raise ValueError(f"`{file}` is not a file.")

    return hmmer_reader.read(file)


def create_profile(hmmfile: hmmer_reader.HMMEReader):
    from ._alphabet import Alphabet
    from ._hmm import HMM
    from ._state import MuteState, NormalState

    alphabet = Alphabet(hmmfile.alphabet)

    hmm = HMM(alphabet)

    D0 = MuteState("D0", alphabet)
    prev = {
        "M": MuteState("B", alphabet),
        "I": NormalState("I0", alphabet, hmmfile.insert(0, True)),
        "D": D0,
    }
    hmm.add_state(prev["M"], LOG(1.0))
    hmm.add_state(prev["I"])
    hmm.add_state(prev["D"])

    for m in range(1, hmmfile.M + 1):
        curr = {
            "M": NormalState(f"M{m}", alphabet, hmmfile.match(m, True)),
            "I": NormalState(f"I{m}", alphabet, hmmfile.insert(m, True)),
            "D": MuteState(f"D{m}", alphabet),
        }
        hmm.add_state(curr["M"])
        hmm.add_state(curr["I"])
        hmm.add_state(curr["D"])

        trans = hmmfile.trans(m - 1, True)

        hmm.set_trans(prev["M"], curr["M"], trans["MM"])
        hmm.set_trans(prev["M"], prev["I"], trans["MI"])
        hmm.set_trans(prev["M"], curr["D"], trans["MD"])
        hmm.set_trans(prev["I"], curr["M"], trans["IM"])
        hmm.set_trans(prev["I"], prev["I"], trans["II"])
        hmm.set_trans(prev["D"], curr["M"], trans["DM"])
        hmm.set_trans(prev["D"], curr["D"], trans["DD"])

        prev = curr

    E = MuteState("E", alphabet)
    hmm.add_state(E)
    hmm.set_trans(E, E, LOG(1.0))

    M = hmmfile.M
    trans = hmmfile.trans(M, True)

    hmm.set_trans(curr["M"], E, trans["MM"])
    hmm.set_trans(curr["M"], curr["I"], trans["MI"])
    hmm.set_trans(curr["I"], E, trans["IM"])
    hmm.set_trans(curr["I"], curr["I"], trans["II"])
    hmm.set_trans(curr["D"], E, trans["DM"])

    hmm.del_state(D0)
    # hmm.rename_state("M0", "B")
    # hmm.normalize()

    return hmm
