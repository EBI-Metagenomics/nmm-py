from typing import List, Dict

GENCODE = {
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
    def __init__(self, name: str = "standard"):
        self._gencode = GENCODE[name]

    def codons(self, amino_acid: bytes) -> List[bytes]:
        amino_acid = amino_acid.upper()
        return self._gencode.get(amino_acid, [])


def generate_codon_lprobs(aa_lprobs: Dict[bytes, float], gencode: GeneticCode):
    from math import log

    codon_lprobs: Dict[bytes, float] = {}
    for aa, logp in aa_lprobs.items():
        codons = gencode.codons(aa)
        if len(codons) == 0:
            continue
        logp_norm = log(len(codons))
        for codon in codons:
            codon_lprobs[codon] = logp - logp_norm

    return codon_lprobs
