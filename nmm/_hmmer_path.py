from typing import Union, List
from ._path import Path

from ._ffi import ffi, lib


class Fragment:
    def __init__(self, seq: bytes, start: int, end: int):
        self._seq = seq
        self._start = start
        self._end = end

    @property
    def sequence(self) -> str:
        return self._seq[self._start : self._end].decode()


class HomoFragment(Fragment):
    def __init__(self, seq: bytes, start: int, end: int):
        super().__init__(seq, start, end)

    @property
    def homologous(self):
        return True

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.sequence}>"


class NonHomoFragment(Fragment):
    def __init__(self, seq: bytes, start: int, end: int):
        super().__init__(seq, start, end)

    @property
    def homologous(self):
        return False

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.sequence}>"


class HMMERResult:
    def __init__(self, score: float, seq: bytes, path: Path):
        self._score = score

        self._fragments: List[Union[HomoFragment, NonHomoFragment]] = []

        frag_start = frag_end = 0
        homologous = False

        step = lib.imm_path_first(path.cdata)
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
            step = lib.imm_path_next(path.cdata, step)

    @property
    def fragments(self):
        return self._fragments

    @property
    def score(self) -> float:
        return self._score
