from typing import Dict, List, Set

from ._nmm import Base, Codon

GENCODE: Dict[str, Dict[bytes, List[bytes]]] = {
    "standard": {
        b"F": [b"UUU", b"UUC"],
        b"L": [b"UUA", b"UUG", b"CUU", b"CUC", b"CUA", b"CUG"],
        b"I": [b"AUU", b"AUC", b"AUA"],
        b"M": [b"AUG"],
        b"V": [b"GUU", b"GUC", b"GUA", b"GUG"],
        b"S": [b"UCU", b"UCC", b"UCA", b"UCG", b"AGU", b"AGC"],
        b"P": [b"CCU", b"CCC", b"CCA", b"CCG"],
        b"T": [b"ACU", b"ACC", b"ACA", b"ACG"],
        b"A": [b"GCU", b"GCC", b"GCA", b"GCG"],
        b"Y": [b"UAU", b"UAC"],
        b"*": [b"UAA", b"UAG", b"UGA"],
        b"H": [b"CAU", b"CAC"],
        b"Q": [b"CAA", b"CAG"],
        b"N": [b"AAU", b"AAC"],
        b"K": [b"AAA", b"AAG"],
        b"D": [b"GAU", b"GAC"],
        b"E": [b"GAA", b"GAG"],
        b"C": [b"UGU", b"UGC"],
        b"W": [b"UGG"],
        b"R": [b"CGU", b"CGC", b"CGA", b"CGG", b"AGA", b"AGG"],
        b"G": [b"GGU", b"GGC", b"GGA", b"GGG"],
    }
}


class GeneticCode:
    def __init__(self, base: Base, name: str = "standard"):

        self._gencode: Dict[bytes, List[Codon]] = {}
        for aa, triplets in GENCODE[name].items():
            gcode = self._gencode[aa]
            for triplet in triplets:
                gcode.append(Codon(triplet, base))

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
