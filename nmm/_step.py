from ._ffi import ffi, lib


class CStep:
    """
    Wrapper around the C implementation of a path step.

    A step is composed of a state and an emitted sequence length. The user should not need to
    directly call the constructor of this class but instead use the methods from the `Path` class.

    Parameters
    ----------
    imm_step : `<cdata 'struct imm_step *'>`
        Step pointer.
    """

    def __init__(self, imm_step: ffi.CData):
        if imm_step == ffi.NULL:
            raise RuntimeError("`imm_step` is NULL.")
        self._imm_step = imm_step

    @property
    def imm_step(self) -> ffi.CData:
        return self._imm_step

    @property
    def seq_len(self) -> int:
        return lib.imm_step_seq_len(self.imm_step)

    def __del__(self):
        if self._imm_step != ffi.NULL:
            lib.imm_step_destroy(self._imm_step)

    def __str__(self) -> str:
        state = lib.imm_step_state(self._imm_step)
        name: str = ffi.string(lib.imm_state_get_name(state))
        return f"{name},{self.seq_len}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
