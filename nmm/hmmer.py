from pathlib import Path
from typing import NamedTuple, Union
from ._state import MuteState, NormalState

import hmmer_reader

from ._log import LOG0


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
    curr = _States(
        M=MuteState("B", alphabet),
        I=NormalState("I0", alphabet, hmmfile.insert(0, True)),
        D=D0,
    )
    curr.I.normalize()
    hmm.add_state(curr.M, 0.0, name="B")
    hmm.add_state(curr.I, name="I0")
    hmm.add_state(curr.D, name="D0")

    for m in range(1, hmmfile.M + 1):
        prev = curr

        curr = _States(
            M=NormalState(f"M{m}", alphabet, hmmfile.match(m, True)),
            I=NormalState(f"I{m}", alphabet, hmmfile.insert(m, True)),
            D=MuteState(f"D{m}", alphabet),
        )
        curr.M.normalize()
        curr.I.normalize()
        hmm.add_state(curr.M, name=f"M{m}")
        hmm.add_state(curr.I, name=f"I{m}")
        hmm.add_state(curr.D, name=f"D{m}")

        trans = hmmfile.trans(m - 1, True)

        hmm.set_trans(prev.M, curr.M, trans["MM"])
        hmm.set_trans(prev.M, prev.I, trans["MI"])
        hmm.set_trans(prev.M, curr.D, trans["MD"])
        hmm.set_trans(prev.I, curr.M, trans["IM"])
        hmm.set_trans(prev.I, prev.I, trans["II"])
        hmm.set_trans(prev.D, curr.M, trans["DM"])
        hmm.set_trans(prev.D, curr.D, trans["DD"])

    E = MuteState("E", alphabet)
    hmm.add_state(E, name="E")
    hmm.set_trans(E, E, 0.0)

    M = hmmfile.M
    trans = hmmfile.trans(M, True)

    hmm.set_trans(curr.M, E, trans["MM"])
    hmm.set_trans(curr.M, curr.I, trans["MI"])
    hmm.set_trans(curr.I, E, trans["IM"])
    hmm.set_trans(curr.I, curr.I, trans["II"])
    hmm.set_trans(curr.D, E, trans["DM"])

    hmm.del_state(D0)
    hmm.normalize()
    hmm.set_trans(E, E, LOG0)

    return hmm


def create_core_profile(hmmfile: hmmer_reader.HMMEReader):
    from ._alphabet import Alphabet
    from ._hmm import HMM
    from ._state import MuteState, NormalState

    alphabet = Alphabet(hmmfile.alphabet)

    hmm = HMM(alphabet)

    curr = _States(
        M=MuteState("B", alphabet),
        I=NormalState("I0", alphabet, hmmfile.insert(0, True)),
        D=MuteState("D0", alphabet),
    )
    curr.I.normalize()

    hmm.add_state(curr.M, 0.0, name="B")
    hmm.add_state(curr.I, name="I0")
    hmm.add_state(curr.D, name="D0")

    for m in range(1, hmmfile.M + 1):
        prev = curr

        curr = _States(
            M=NormalState(f"M{m}", alphabet, hmmfile.match(m, True)),
            I=NormalState(f"I{m}", alphabet, hmmfile.insert(m, True)),
            D=MuteState(f"D{m}", alphabet),
        )
        curr.M.normalize()
        curr.I.normalize()
        hmm.add_state(curr.M, name=f"M{m}")
        hmm.add_state(curr.I, name=f"I{m}")
        hmm.add_state(curr.D, name=f"D{m}")

        trans = hmmfile.trans(m - 1, True)

        hmm.set_trans(prev.M, curr.M, trans["MM"])
        hmm.set_trans(prev.M, prev.I, trans["MI"])
        hmm.set_trans(prev.M, curr.D, trans["MD"])
        hmm.set_trans(prev.I, curr.M, trans["IM"])
        hmm.set_trans(prev.I, prev.I, trans["II"])
        hmm.set_trans(prev.D, curr.M, trans["DM"])
        hmm.set_trans(prev.D, curr.D, trans["DD"])

    E = MuteState("E", alphabet)
    hmm.add_state(E, name="E")

    M = hmmfile.M
    trans = hmmfile.trans(M, True)

    hmm.set_trans(curr.M, E, trans["MM"])
    hmm.set_trans(curr.M, curr.I, trans["MI"])
    hmm.set_trans(curr.I, E, trans["IM"])
    hmm.set_trans(curr.I, curr.I, trans["II"])
    hmm.set_trans(curr.D, E, trans["DM"])

    hmm.del_state(hmm.find_state("D0"))

    hmm.set_trans(E, E, 0.0)
    hmm.normalize()
    hmm.set_trans(E, E, LOG0)

    return hmm


def create_global_profile(hmm):
    from ._state import MuteState, NormalState

    X = 1.0

    alphabet = hmm.alphabet
    B = hmm.find_state("B")

    hmm.set_start_lprob(B, LOG0)

    S = MuteState("S", alphabet)
    hmm.add_state(S, 0.0, "S")

    I1: NormalState = hmm.find_state("I1")
    bg_emission_table = I1.emission_table()

    N = NormalState("N", alphabet, bg_emission_table)
    hmm.add_state(N, name="N")
    hmm.set_trans(S, N, 0.0)
    hmm.set_trans(N, N, X)
    hmm.set_trans(N, B, 0.0)

    S = MuteState("S", alphabet)
    pass


_States = NamedTuple(
    "States",
    [("M", Union[NormalState, MuteState]), ("I", NormalState), ("D", MuteState)],
)
