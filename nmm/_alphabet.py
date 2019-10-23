from ._ffi import ffi, lib


class Alphabet:
    def __init__(self, symbols: str):
        self._abc = ffi.NULL
        self._symbols = symbols.encode()
        self._abc = lib.imm_abc_create(self._symbols)
        if self._abc == ffi.NULL:
            raise RuntimeError("Could not create alphabet.")

    @property
    def cdata(self):
        return self._abc

    @property
    def length(self) -> int:
        return lib.imm_abc_length(self._abc)

    def has_symbol(self, symbol_id: str) -> bool:
        return lib.imm_abc_has_symbol(self._abc, symbol_id.encode()) == 1

    def symbol_idx(self, symbol_id: str) -> int:
        return lib.imm_abc_symbol_idx(self._abc, symbol_id.encode())

    def symbol_id(self, symbol_idx: int) -> str:
        return lib.imm_abc_symbol_id(self._abc, symbol_idx).decode()

    @property
    def symbols(self) -> str:
        return self._symbols.decode()

    def __str__(self) -> str:
        return f"{{{self.symbols}}}"

    def __repr__(self) -> str:
        return f"{{{self.__class__.__name__}:{self.symbols}}}"

    def __del__(self):
        if self._abc != ffi.NULL:
            lib.imm_abc_destroy(self._abc)
