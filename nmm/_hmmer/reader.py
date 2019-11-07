import pathlib
from io import TextIOBase
from typing import Union

import hmmer_reader


def read_hmmer(file: Union[str, pathlib.Path, TextIOBase]) -> hmmer_reader.HMMEReader:

    if isinstance(file, str):
        file = pathlib.Path(file)

    if isinstance(file, pathlib.Path):
        if not file.exists():
            raise FileNotFoundError(f"`{file}` does not exist.")

        if not file.is_file():
            raise ValueError(f"`{file}` is not a file.")

    return hmmer_reader.read(file)
