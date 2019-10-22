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

    breakpoint()
    hmm.add_state(MuteState("M0", alphabet), LOG(1.0))
    hmm.add_state(NormalState("I0", alphabet, hmmfile.insert(0, True)))
    hmm.add_state(MuteState("D0", alphabet))

    for m in range(1, hmmfile.M + 1):
        hmm.add_state(NormalState(f"M{m}", hmmfile.match(m, True)))
        hmm.add_state(NormalState(f"I{m}", hmmfile.insert(m, True)))
        hmm.add_state(MuteState(f"D{m}", alphabet, False))

        trans = hmmfile.trans(m - 1, True)

        hmm.set_trans(f"M{m-1}", f"M{m}", trans["MM"])
        hmm.set_trans(f"M{m-1}", f"I{m-1}", trans["MI"])
        hmm.set_trans(f"M{m-1}", f"D{m}", trans["MD"])
        hmm.set_trans(f"I{m-1}", f"M{m}", trans["IM"])
        hmm.set_trans(f"I{m-1}", f"I{m-1}", trans["II"])
        hmm.set_trans(f"D{m-1}", f"M{m}", trans["DM"])
        hmm.set_trans(f"D{m-1}", f"D{m}", trans["DD"])

    hmm.add_state(MuteState("E", alphabet, True))
    hmm.set_trans("E", "E", LOG(1.0))

    M = hmmfile.M
    trans = hmmfile.trans(M, True)
    hmm.set_trans(f"M{M}", f"E", trans["MM"])
    hmm.set_trans(f"M{M}", f"I{M}", trans["MI"])
    hmm.set_trans(f"I{M}", f"E", trans["IM"])
    hmm.set_trans(f"I{M}", f"I{M}", trans["II"])
    hmm.set_trans(f"D{M}", f"E", trans["DM"])

    hmm.delete_state("D0")
    hmm.rename_state("M0", "B")
    hmm.normalize()
    return hmm
