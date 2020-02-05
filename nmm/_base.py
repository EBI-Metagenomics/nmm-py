from ._alphabet import CAlphabet
from ._ffi import ffi, lib


class CBase:
    def __init__(self, nmm_base: ffi.CData):
        super().__init__()
        if nmm_base == ffi.NULL:
            raise RuntimeError("`nmm_base` is NULL.")
        self._nmm_base = nmm_base

    @property
    def nmm_base(self) -> ffi.CData:
        return self._nmm_base

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_abc_symbols(lib.nmm_base_get_abc(self._nmm_base)))

    def __del__(self):
        if self._nmm_base != ffi.NULL:
            lib.nmm_base_destroy(self._nmm_base)

    def __str__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{symbols}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Base(CBase):
    """
    Base is a nucleotide letter.

    base : Union[bytes, str, int]
        A single letter.
    """

    def __init__(self, calphabet: CAlphabet):
        self._calphabet = calphabet

        nmm_base = lib.nmm_base_create(calphabet.imm_abc)
        if nmm_base == ffi.NULL:
            raise RuntimeError("`nmm_base_create` failed.")

        super().__init__(nmm_base)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
