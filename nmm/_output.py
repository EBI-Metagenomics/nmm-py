from __future__ import annotations

from typing import Type

from ._cdata import CData
from ._ffi import ffi, lib
from ._model import Model

__all__ = ["Output"]


class Output:
    def __init__(self, nmm_output: CData):
        if nmm_output == ffi.NULL:
            raise RuntimeError("`nmm_output` is NULL.")
        self._nmm_output = nmm_output

    @classmethod
    def create(cls: Type[Output], filepath: bytes) -> Output:
        return cls(lib.nmm_output_create(filepath))

    def write(self, model: Model):
        err: int = lib.nmm_output_write(self._nmm_output, model.nmm_model)
        if err != 0:
            raise RuntimeError("Could not write model.")

    def close(self):
        err: int = lib.nmm_output_close(self._nmm_output)
        if err != 0:
            raise RuntimeError("Could not close output.")

    def __del__(self):
        if self._nmm_output != ffi.NULL:
            self.close()
            lib.nmm_output_destroy(self._nmm_output)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()
