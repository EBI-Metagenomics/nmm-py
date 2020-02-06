from ._ffi import ffi, lib


class CAlphabet:
    """
    Wrapper around the C implementation of alphabet.

    Parameters
    ----------
    imm_abc : `<cdata 'struct imm_abc *'>`.
    """

    def __init__(self, imm_abc: ffi.CData):
        super().__init__()
        if imm_abc == ffi.NULL:
            raise RuntimeError("`imm_abc` is NULL.")
        self._imm_abc = imm_abc

    @property
    def imm_abc(self) -> ffi.CData:
        return self._imm_abc

    @property
    def length(self) -> int:
        return lib.imm_abc_length(self._imm_abc)

    @property
    def symbols(self) -> bytes:
        return ffi.string(lib.imm_abc_symbols(self._imm_abc))

    def has_symbol(self, symbol_id: bytes) -> bool:
        return lib.imm_abc_has_symbol(self._imm_abc, symbol_id)

    def symbol_idx(self, symbol_id: bytes) -> int:
        return lib.imm_abc_symbol_idx(self._imm_abc, symbol_id)

    def symbol_id(self, symbol_idx: int) -> bytes:
        return lib.imm_abc_symbol_id(self._imm_abc, symbol_idx)

    def __del__(self):
        if self._imm_abc != ffi.NULL:
            lib.imm_abc_destroy(self._imm_abc)

    def __str__(self) -> str:
        return f"{{{self.symbols.decode()}}}"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"


class Alphabet(CAlphabet):
    """
    Alphabet set for Markov models.

    Parameters
    ----------
    symbols : bytes
        Set of symbols as an array of bytes.
    """

    def __init__(self, symbols: bytes, any_symbol: bytes):
        if len(any_symbol) != 1:
            raise ValueError("`any_symbol` has length different than 1.")

        imm_abc = lib.imm_abc_create(symbols, any_symbol)
        if imm_abc == ffi.NULL:
            raise RuntimeError("`imm_abc_create` failed.")

        super().__init__(imm_abc)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}:{str(self)}>"
