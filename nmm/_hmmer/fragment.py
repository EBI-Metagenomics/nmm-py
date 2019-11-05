class Fragment:
    def __init__(self, seq: bytes, start: int, end: int):
        self._seq = seq
        self._start = start
        self._end = end

    @property
    def sequence(self) -> bytes:
        return self._seq[self._start : self._end]


class HomoFragment(Fragment):
    def __init__(self, seq: bytes, start: int, end: int):
        super().__init__(seq, start, end)

    @property
    def homologous(self):
        return True

    def __repr__(self):
        seq = self.sequence.encode()
        return f"<{self.__class__.__name__}:{seq}>"


class NonHomoFragment(Fragment):
    def __init__(self, seq: bytes, start: int, end: int):
        super().__init__(seq, start, end)

    @property
    def homologous(self):
        return False

    def __repr__(self):
        seq = self.sequence.encode()
        return f"<{self.__class__.__name__}:{seq}>"
