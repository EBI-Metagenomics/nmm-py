from typing import List, Union

from .._ffi import ffi, lib
from .._path import CPath
from .fragment import HomoFragment, NonHomoFragment


class Result:
    def __init__(self, score: float, seq: bytes, path: CPath):
        self._score = score

        self._fragments: List[Union[HomoFragment, NonHomoFragment]] = []

        frag_start = frag_end = 0
        homologous = False

        step = lib.imm_path_first(path.imm_path)
        while step != ffi.NULL:
            cname = lib.imm_state_get_name(lib.imm_step_state(step))
            name = ffi.string(cname).decode()
            seq_len = lib.imm_step_seq_len(step)

            if not homologous and name.startswith("M"):
                if frag_start < frag_end:
                    self._fragments.append(NonHomoFragment(seq, frag_start, frag_end))
                homologous = True
                frag_start = frag_end

            elif homologous and name.startswith("E"):
                if frag_start < frag_end:
                    self._fragments.append(HomoFragment(seq, frag_start, frag_end))
                homologous = False
                frag_start = frag_end

            frag_end += seq_len
            step = lib.imm_path_next(path.imm_path, step)

    @property
    def fragments(self):
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
