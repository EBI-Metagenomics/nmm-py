from __future__ import annotations

from typing import Type

from .._cdata import CData
from .._ffi import ffi, lib
from ._model import Model


class Output:
    def __init__(self, nmm_output: CData):
        if nmm_output == ffi.NULL:
            raise RuntimeError("`nmm_output` is NULL.")
        self._nmm_output = nmm_output

    @classmethod
    def create(cls: Type[Output], filepath: bytes) -> Output:
        return cls(lib.nmm_output_create(filepath))

    def write(self, model: Model):
        lib.nmm_output_write(self._nmm_output, model.nmm_model)

    def __del__(self):
        if self._nmm_output != ffi.NULL:
            lib.nmm_output_destroy(self._nmm_output)
