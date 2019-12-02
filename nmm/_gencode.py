from typing import List, Dict, Set
from ._codon import Codon

GENCODE: Dict[str, Dict[bytes, List[Codon]]] = {
    "standard": {
        b"F": [Codon("UUU"), Codon("UUC")],
        b"L": [
            Codon("UUA"),
            Codon("UUG"),
            Codon("CUU"),
            Codon("CUC"),
            Codon("CUA"),
            Codon("CUG"),
        ],
        b"I": [Codon("AUU"), Codon("AUC"), Codon("AUA")],
        b"M": [Codon("AUG")],
        b"V": [Codon("GUU"), Codon("GUC"), Codon("GUA"), Codon("GUG")],
        b"S": [
            Codon("UCU"),
            Codon("UCC"),
            Codon("UCA"),
            Codon("UCG"),
            Codon("AGU"),
            Codon("AGC"),
        ],
        b"P": [Codon("CCU"), Codon("CCC"), Codon("CCA"), Codon("CCG")],
        b"T": [Codon("ACU"), Codon("ACC"), Codon("ACA"), Codon("ACG")],
        b"A": [Codon("GCU"), Codon("GCC"), Codon("GCA"), Codon("GCG")],
        b"Y": [Codon("UAU"), Codon("UAC")],
        b"*": [Codon("UAA"), Codon("UAG"), Codon("UGA")],
        b"H": [Codon("CAU"), Codon("CAC")],
        b"Q": [Codon("CAA"), Codon("CAG")],
        b"N": [Codon("AAU"), Codon("AAC")],
        b"K": [Codon("AAA"), Codon("AAG")],
        b"D": [Codon("GAU"), Codon("GAC")],
        b"E": [Codon("GAA"), Codon("GAG")],
        b"C": [Codon("UGU"), Codon("UGC")],
        b"W": [Codon("UGG")],
        b"R": [
            Codon("CGU"),
            Codon("CGC"),
            Codon("CGA"),
            Codon("CGG"),
            Codon("AGA"),
            Codon("AGG"),
        ],
        b"G": [Codon("GGU"), Codon("GGC"), Codon("GGA"), Codon("GGG")],
    }
}


class GeneticCode:
    def __init__(self, name: str = "standard"):
        self._gencode = GENCODE[name]
        self._amino_acid: Dict[Codon, bytes] = {}
        for aa, codons in self._gencode.items():
            for codon in codons:
                self._amino_acid[codon] = aa

    def codons(self, amino_acid: bytes) -> List[Codon]:
        amino_acid = amino_acid.upper()
        return self._gencode.get(amino_acid, [])

    def amino_acid(self, codon: Codon) -> bytes:
        return self._amino_acid[codon]

    def amino_acids(self) -> Set[bytes]:
        return set(self._gencode.keys())
