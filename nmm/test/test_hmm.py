from math import log

from imm import HMM, MuteState, Sequence
from imm.testing import assert_allclose

from nmm import BaseAlphabet, BaseTable, Codon, CodonProb, CodonTable, FrameState


def test_hmm():
    abc = BaseAlphabet.create(b"ACGU", b"X")
    baset = BaseTable.create(abc, (log(0.25), log(0.25), log(0.25), log(0.25)))

    codonp = CodonProb.create(abc)
    codonp.set_lprob(Codon.create(b"AUG", abc), log(0.8))
    codonp.set_lprob(Codon.create(b"AUU", abc), log(0.1))

    B = MuteState.create(b"B", abc)
    M1 = FrameState.create(b"M1", baset, CodonTable.create(codonp), 0.02)
    M2 = FrameState.create(b"M2", baset, CodonTable.create(codonp), 0.01)
    E = MuteState.create(b"E", abc)

    hmm = HMM.create(abc)
    hmm.add_state(B, log(0.5))
    hmm.add_state(M1)
    hmm.add_state(M2)
    hmm.add_state(E)

    hmm.set_transition(B, M1, log(0.8))
    hmm.set_transition(B, M2, log(0.2))
    hmm.set_transition(M1, M2, log(0.1))
    hmm.set_transition(M1, E, log(0.4))
    hmm.set_transition(M2, E, log(0.3))

    dp = hmm.create_dp(E)
    results = dp.viterbi(Sequence.create(b"AUGAUU", abc))
    assert len(results) == 1
    assert_allclose(results[0].loglikelihood, -7.069201008427531)
