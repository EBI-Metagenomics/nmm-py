from io import TextIOBase


class FASTAParserError(Exception):
    pass


class FASTAReader:
    def __init__(self, file: TextIOBase):
        self._file = file

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._file.close()
        del type
        del value
        del traceback

    def _next_defline(self):
        while True:
            line = self._file.readline()
            if line == "":
                return None

            line = line.strip()
            if line.startswith(">"):
                return line
            elif line != "":
                raise FASTAParserError()

    def _next_sequence(self):
        while True:
            line = self._file.readline()
            if line == "":
                raise FASTAParserError()

            line = line.strip()
            if not line.startswith(">"):
                return line
            elif line != "":
                raise FASTAParserError()

    def items(self):

        while True:
            defline = self._next_defline()
            if defline is None:
                return

            sequence = self._next_sequence()
            yield (defline, sequence)

