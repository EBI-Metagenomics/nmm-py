from abc import ABC, abstractmethod

from ._ffi import ffi, lib
from ._state import CState, StateBase


class StepBase(ABC):
    """
    Path step.

    A step is composed of a state and an emitted sequence length. The user should not need to
    directly call the constructor of this class but instead use the methods from the `Path` class.
    """

    @property
    @abstractmethod
    def state(self) -> StateBase:
        raise NotImplementedError()

    @property
    @abstractmethod
    def seq_len(self) -> int:
        raise NotImplementedError()

    def __str__(self) -> str:
        name = self.state.name.decode()
        return f"<{name},{self.seq_len}>"

    def __repr__(self) -> str:
        name = self.state.name.decode()
        return f"<{self.__class__.__name__}:{name},{self.seq_len}>"


class CStep(StepBase):
    """
    Wrapper around the C implementation of path step.

    Parameters
    ----------
    imm_step : ffi.CData
        A non-null `<cdata 'struct imm_step *'>`.
    """

    def __init__(self, imm_step: ffi.CData):
        super().__init__()
        if imm_step == ffi.NULL:
            raise RuntimeError("`imm_step` is NULL.")
        self._imm_step = imm_step

    @property
    def imm_step(self) -> ffi.CData:
        return self._imm_step

    @property
    def state(self) -> CState:
        return CState(lib.imm_step_state(self.imm_step))

    @property
    def seq_len(self) -> int:
        return lib.imm_step_seq_len(self.imm_step)

    def __repr__(self) -> str:
        name = self.state.name.decode()
        return f"<{self.__class__.__name__}:{name},{self.seq_len}>"
