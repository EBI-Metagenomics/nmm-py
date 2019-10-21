from math import exp, isinf, inf, isnan
from ._log import LOG
from ._state import State
from ._alphabet import Alphabet
import functools

from ._ffi import ffi, lib


class HMM:
    def __init__(self, alphabet: Alphabet):
        self._alphabet = alphabet
        self._hmm = lib.imm_hmm_create(self._alphabet.cdata)
        self._states = {}

    def __del__(self):
        if self._hmm != ffi.NULL:
            lib.imm_hmm_destroy(self._hmm)

    @property
    def states(self):
        return self._states

    def set_start_lprob(self, state: State, lprob: float):
        err: int = lib.imm_hmm_set_start_lprob(self._hmm, state.cdata, lprob)
        if err != 0:
            raise ValueError("Could not set start probability.")

    def trans(self, a: State, b: State):
        """
        Parameters
        ----------
        a : State
            Source state.
        b : State
            Destination state.
        """
        lprob: float = lib.imm_hmm_get_trans(self._hmm, a.cdata, b.cdata)
        if isnan(lprob):
            raise ValueError("Could not retrieve transition probability.")
        return lprob

    def set_trans(self, a: State, b: State, lprob: float):
        """
        Parameters
        ----------
        a : State
            Source state name.
        b : State
            Destination state name.
        lprob : float
            Transition probability in log-space.
        """
        err: int = lib.imm_hmm_set_trans(self._hmm, a.cdata, b.cdata, lprob)
        if err != 0:
            raise ValueError("Could not set transition probability.")

    @property
    def alphabet(self):
        return self._alphabet

    def add_state(self, state: State, start_lprob: float = LOG(0.0)):
        """
        Parameters
        ----------
        state
            Add state.
        start_lprob : bool
            Log-space probability of being the initial state.
        """
        err: int = lib.imm_hmm_add_state(self._hmm, state.cdata, start_lprob)
        if err != 0:
            raise ValueError("Could not add state %s.", state)
        self._states[state.cdata] = state

    def rename_state(self, old_name: str, new_name: str):
        if old_name not in self._states:
            raise ValueError(f"State name `{old_name}` does not exist.")

        if new_name in self._states:
            raise ValueError(f"State name `{new_name}` already exists.")

        self._states[new_name] = self._states.pop(old_name)
        self._states[new_name].name = new_name
        self._init_logps[new_name] = self._init_logps.pop(old_name)

        for k, v in self._trans.items():
            if old_name in v:
                v[new_name] = v.pop(old_name)
        self._trans[new_name] = self._trans.pop(old_name)

    def delete_state(self, name):
        if name not in self._states:
            raise ValueError(f"State name `{name}` does not exist.")

        del self._states[name]
        del self._init_logps[name]
        del self._trans[name]
        for v in self._trans.values():
            if name in v:
                del v[name]

    def normalize(self):
        err: int = lib.imm_hmm_normalize(self._hmm)
        if err != 0:
            raise ValueError("Normalization error.")

    def likelihood(self, seq: str, state_path: list, log_space: bool = False):
        if len(state_path) == 0:
            if len(seq) == 0:
                if log_space:
                    return LOG(1.0)
                return 1.0
            if log_space:
                return LOG(0.0)
            return 0.0

        self._assure_states_exist([i[0] for i in state_path])
        head = state_path[0]
        qt = self._states[head[0]]
        ft = head[1]
        if ft > len(seq):
            return 0.0
        logp = self.init_prob(qt.name, True) + qt.prob(seq[:ft], True)

        seq = seq[ft:]
        qt_1 = qt
        for head in state_path[1:]:
            qt = self._states[head[0]]
            ft = head[1]
            logp += qt.prob(seq[:ft], True) + self.trans(qt_1.name, qt.name, True)
            seq = seq[ft:]
            qt_1 = qt

        if log_space:
            return logp
        return exp(logp)

    def draw(self, filepath, emissions=0, init_prob=True, digits=3, view=False):
        from graphviz import Digraph

        graph = Digraph()

        for state in self._states.values():
            shape = "circle"

            if init_prob:
                p = self.init_prob(state.name, log_space=False)
                p = round(p, digits)
                if p > 0:
                    state_label = f"{state.name}: {p}"
                else:
                    state_label = f"{state.name}"
            else:
                state_label = f"{state.name}"

            if emissions > 0:
                emission = state.emission(log_space=False)
                emission = emission[:emissions]
                label = _format_emission_table(emission, state_label, digits)
            else:
                label = state_label

            graph.node(state.name, label, shape=shape)

        for state0, trans in self._trans.items():
            for state1, logp in trans.items():
                p = exp(logp)
                p = round(p, digits)
                if p > 0:
                    graph.edge(state0, state1, label=f"{p}")

        graph.render(filepath, view=view)

    def viterbi(self, seq: str, end_state: str, log_space: bool = False):
        self._check_mute_cycle()
        self._vcache = {}
        max_logp = LOG(0.0)
        best_path = []

        end_state = self._states[end_state]
        for ft in range(end_state.min_len, end_state.max_len + 1):
            tup = self._viterbi(seq, end_state, ft)
            if tup[0] > max_logp:
                max_logp = tup[0]
                best_path = tup[1] + [(end_state, ft)]

        best_path = [(qt.name, ft) for qt, ft in best_path]
        if log_space:
            return max_logp, best_path
        return exp(max_logp), best_path

    @functools.lru_cache(maxsize=65536)
    def _viterbi(self, seq: str, qt: str, ft: int):
        max_logp = LOG(0.0)
        best_path = []
        if ft > len(seq):
            return max_logp, best_path

        emission_prob = qt.prob(seq[len(seq) - ft :], True)
        if emission_prob <= max_logp:
            return max_logp, best_path

        seq_end = len(seq) - ft
        prefix = seq[:seq_end]
        for qt_1 in self._states.values():
            if qt_1.min_len > len(prefix):
                continue

            T = self.trans(qt_1.name, qt.name, True)
            if T + emission_prob <= max_logp:
                continue

            for ft_1 in range(qt_1.min_len, qt_1.max_len + 1):
                tup = self._viterbi(prefix, qt_1, ft_1)
                tup = (tup[0] + T + emission_prob, tup[1] + [(qt_1, ft_1)])

                if tup[0] > max_logp:
                    max_logp = tup[0]
                    best_path = tup[1]

        if len(seq) - ft == 0:
            v = emission_prob + self.init_prob(qt.name, True)
            if v > max_logp:
                max_logp = v
                best_path = []

        return max_logp, best_path

    def _draw_initial_state(self, random):
        names = self._init_logps.keys()
        probs = [exp(v) for v in self._init_logps.values()]
        name = random.choice(list(names), p=probs)
        return self._states[name]

    def _transition(self, state: State, random):
        trans = self._trans[state.name]
        names = trans.keys()
        probs = [exp(v) for v in trans.values()]
        name = random.choice(list(names), p=probs)
        return self._states[name]

    def _normalize_trans(self):
        from scipy.special import logsumexp

        names = self._states.keys()
        nstates = len(names)
        for a in names:
            logprobs = list(self._trans[a].values())

            if len(logprobs) == 0:
                continue
                # for b in names:
                #     self._trans[a][b] = -LOG(nstates)
            # else:
            logprob_norm = logsumexp(logprobs)
            if isinf(logprob_norm):
                for b in names:
                    self._trans[a][b] = -LOG(nstates)
            else:
                for b in self._trans[a].keys():
                    self._trans[a][b] -= logprob_norm

    def _normalize_init_logps(self):
        from scipy.special import logsumexp

        logprobs = list(self._init_logps.values())
        names = self._states.keys()
        nstates = len(names)

        logp_norm = logsumexp(logprobs)
        if isinf(logp_norm):
            for a in names:
                self._init_logps[a] = -LOG(nstates)
        else:
            for a in self._init_logps.keys():
                self._init_logps[a] -= logp_norm

    def _assure_states_exist(self, states):
        for state in states:
            if state not in self._states:
                raise ValueError(f"State `{state}` does not exist.")

    def _check_mute_cycle(self):
        state_marks = {name: "initial" for name in self._states}
        for state_name in self._states:
            self._mute_cycle_visit(state_name, state_marks)

    def _mute_cycle_visit(self, state: str, state_marks):
        if self._states[state].min_len > 0:
            return
        if state_marks[state] == "permament":
            return
        if state_marks[state] == "temporary":
            raise ValueError("Mute cycles are not allowed.")

        state_marks[state] = "temporary"
        for dst in self._trans[state]:

            if self._trans[state][dst] > float("-inf"):
                self._mute_cycle_visit(dst, state_marks)

        state_marks[state] = "permanent"


def _format_emission_table(emission, name, digits):
    rows = ""
    for row in emission:
        seq = row[0]
        p = round(row[1], digits)
        if p == 0.0:
            break
        rows += f"<TR><TD>{seq}</TD><TD>{p}</TD></TR>"

    tbl_fmt = "BORDER='0' CELLBORDER='1' "
    tbl_fmt += "CELLSPACING='0' CELLPADDING='4'"
    tbl_str = f"<<TABLE {tbl_fmt}>"
    tbl_str += f"<TR><TD COLSPAN='2'>{name}</TD></TR>"
    tbl_str += rows
    tbl_str += "</TABLE>>"

    return tbl_str
