from ._base import Base
from ._string import make_sure_bytes
from ._codon import Codon
from ._state import State

from ._ffi import ffi, lib


class FrameState(State):
    def __init__(self, name: str, base: Base, codon: Codon, epsilon: float):
        """
        Parameters
        ----------
        name : str
            Name.
        """
        super(FrameState, self).__init__(codon.alphabet)
        self._base = base
        self._codon = codon

        name = make_sure_bytes(name)
        self._state = lib.nmm_frame_state_create(name, base.cdata, codon.cdata, epsilon)

    def __del__(self):
        if self._state != ffi.NULL:
            lib.nmm_frame_state_destroy(self._state)

    def __repr__(self):
        return f"<{self.__class__.__name__}:{self.name}>"
