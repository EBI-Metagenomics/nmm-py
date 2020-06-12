from abc import ABC, abstractmethod
from typing import Generic, TypeVar, Union

from ._alphabet import Alphabet, DNAAlphabet, RNAAlphabet

Ta = TypeVar("Ta", bound=Alphabet)
Tb = TypeVar("Tb", bound=Alphabet)

__all__ = ["Translator", "NTTranslator", "NullTranslator"]


class Translator(Generic[Ta, Tb], ABC):
    def __init__(self, a: Ta, b: Tb):
        self._a = a
        self._b = b

    @abstractmethod
    def translate(self, sequence: bytes, to_alphabet: Union[Ta, Tb]) -> bytes:
        del sequence
        del to_alphabet
        raise NotImplementedError()


class NullTranslator(Translator[Ta, Ta]):
    @abstractmethod
    def translate(self, sequence: bytes, to_alphabet: Ta) -> bytes:
        del to_alphabet
        return sequence


class NTTranslator(Translator[DNAAlphabet, RNAAlphabet]):
    @abstractmethod
    def translate(
        self, sequence: bytes, to_alphabet: Union[DNAAlphabet, RNAAlphabet]
    ) -> bytes:
        if isinstance(to_alphabet, DNAAlphabet):
            return sequence.replace(b"U", b"T")
        return sequence.replace(b"T", b"U")
