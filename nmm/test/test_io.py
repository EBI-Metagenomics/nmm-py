from math import log

from nmm import HMM
from nmm.alphabet import Alphabet
from nmm.io import Model, Output
from nmm.prob import lprob_zero
from nmm.state import MuteState, NormalState


def test_nmm_model():

    alphabet = Alphabet.create(b"AC", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, lprob_zero())

    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2)])
    hmm.add_state(M1, lprob_zero())

    M2 = NormalState(b"M2", alphabet, [log(0.4), log(0.6)])
    hmm.add_state(M2, lprob_zero())

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()
    hmm.set_transition(E, E, lprob_zero())

    dp = hmm.create_dp(E)

    model = Model.create(hmm, dp)
    output = Output.create(b"model1.nmm")
    output.write(model)
