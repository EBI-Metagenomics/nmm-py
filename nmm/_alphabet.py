from ._ffi import ffi, lib


class Alphabet:
    def __init__(self, symbols: bytes):
        self._abc = ffi.NULL
        self._symbols = symbols
        self._abc = lib.imm_abc_create(self._symbols)
        if self._abc == ffi.NULL:
            raise RuntimeError("Could not create alphabet.")

    @property
    def cdata(self):
        return self._abc

    @property
    def length(self) -> int:
        return lib.imm_abc_length(self._abc)

    def has_symbol(self, symbol_id: bytes) -> bool:
        return lib.imm_abc_has_symbol(self._abc, symbol_id) == 1

    def symbol_idx(self, symbol_id: bytes) -> int:
        return lib.imm_abc_symbol_idx(self._abc, symbol_id)

    def symbol_id(self, symbol_idx: int) -> bytes:
        return lib.imm_abc_symbol_id(self._abc, symbol_idx)

    @property
    def symbols(self) -> bytes:
        return self._symbols

    def __str__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{symbols}}}"

    def __repr__(self) -> str:
        symbols = self.symbols.decode()
        return f"{{{self.__class__.__name__}:{symbols}}}"

    def __del__(self):
        if self._abc != ffi.NULL:
            lib.imm_abc_destroy(self._abc)
