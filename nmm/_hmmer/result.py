from typing import List, Union

from .._ffi import ffi, lib
from .._path import CPath
from .._step import CStep
from .fragment import HomoFragment, NonHomoFragment


class Result:
    def __init__(self, score: float, seq: bytes, path: CPath):
        self._score = score

        self._fragments: List[Union[HomoFragment, NonHomoFragment]] = []

        frag_start = frag_end = 0
        idx_start = idx_end = 0
        homologous = False

        steps = list(path.steps())
        for step in steps:
            name = step.state.name
            seq_len = step.seq_len

            if not homologous and name.startswith(b"M"):
                if frag_start < frag_end:
                    s = seq[frag_start:frag_end]
                    self._fragments.append(NonHomoFragment(s, steps[idx_start:idx_end]))
                homologous = True
                frag_start = frag_end
                idx_start = idx_end

            elif homologous and name.startswith(b"E"):
                if frag_start < frag_end:
                    s = seq[frag_start:frag_end]
                    self._fragments.append(HomoFragment(s, steps[idx_start:idx_end]))
                homologous = False
                frag_start = frag_end
                idx_start = idx_end

            # subpath.append(step.state.imm_state, step.seq_len)
            frag_end += seq_len
            idx_end += 1
            # step = lib.imm_path_next(path.imm_path, step)

        # step = lib.imm_path_first(path.imm_path)
        # while step != ffi.NULL:
        #     cname = lib.imm_state_get_name(lib.imm_step_state(step))
        #     name = ffi.string(cname)
        #     seq_len = lib.imm_step_seq_len(step)

        #     if not homologous and name.startswith(b"M"):
        #         if frag_start < frag_end:
        #             self._fragments.append(NonHomoFragment(seq, frag_start, frag_end))
        #         homologous = True
        #         frag_start = frag_end

        #     elif homologous and name.startswith(b"E"):
        #         if frag_start < frag_end:
        #             self._fragments.append(HomoFragment(seq, frag_start, frag_end))
        #         homologous = False
        #         frag_start = frag_end

        #     frag_end += seq_len
        #     step = lib.imm_path_next(path.imm_path, step)

    @property
    def fragments(self):
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
