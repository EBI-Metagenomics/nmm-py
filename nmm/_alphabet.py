from ._ffi import ffi, lib


class CAlphabet:
    """
    Wrapper around the C implementation of alphabet set.

    Parameters
    ----------
    cdata : `<cdata 'struct imm_abc *'>`.
    """

    def __init__(self, cdata: ffi.CData):
        self.__cdata = cdata

    @property
    def imm_abc(self) -> ffi.CData:
        return self.__cdata

    @property
    def length(self) -> int:
        return lib.imm_abc_length(self.__cdata)

    def has_symbol(self, symbol_id: bytes) -> bool:
        return lib.imm_abc_has_symbol(self.__cdata, symbol_id) == 1

    def symbol_idx(self, symbol_id: bytes) -> int:
        return lib.imm_abc_symbol_idx(self.__cdata, symbol_id)

    def symbol_id(self, symbol_idx: int) -> bytes:
        return lib.imm_abc_symbol_id(self.__cdata, symbol_idx)

    def __del__(self):
        if self.__cdata != ffi.NULL:
            lib.imm_abc_destroy(self.__cdata)


class Alphabet(CAlphabet):
    """
    Alphabet set for Markov models.

    Parameters
    ----------
    symbols : set of symbols as an array of bytes.
    """

    def __init__(self, symbols: bytes):
        self._symbols = symbols
        cdata = ffi.NULL
        cdata = lib.imm_abc_create(self._symbols)
        if cdata == ffi.NULL:
            raise RuntimeError("`imm_abc_create` failed.")
        super().__init__(cdata)

    @property
    def symbols(self) -> bytes:
        return self._symbols

    def __str__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{symbols}}}"

    def __repr__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{self.__class__.__name__}:{symbols}}}"
