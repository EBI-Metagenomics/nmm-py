from ._string import make_sure_bytes
from ._ffi import ffi, lib


class Alphabet:
    def __init__(self, symbols: bytes):
        self._abc = ffi.NULL
        self._abc = lib.imm_abc_create(make_sure_bytes(symbols))
        if self._abc == ffi.NULL:
            raise RuntimeError(f"Could not create alphabet.")

    def __del__(self):
        if self._abc != ffi.NULL:
            lib.imm_abc_destroy(self._abc)

    @property
    def length(self):
        return lib.imm_abc_length(self._abc)

    def has_symbol(self, symbol_id: bytes):
        return lib.imm_abc_has_symbol(self._abc, make_sure_bytes(symbol_id))

    def symbol_idx(self, symbol_id: bytes):
        return lib.imm_abc_symbol_idx(self._abc, make_sure_bytes(symbol_id))

    def symbol_id(self, symbol_idx: int):
        return lib.imm_abc_symbol_id(self._abc, symbol_idx).decode()

