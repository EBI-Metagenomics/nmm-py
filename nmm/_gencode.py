from typing import List, Dict

GENCODE = {
    "standard": {
        "F": ["UUU", "UUC"],
        "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "I": ["AUU", "AUC", "AUA"],
        "M": ["AUG"],
        "V": ["GUU", "GUC", "GUA", "GUG"],
        "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "P": ["CCU", "CCC", "CCA", "CCG"],
        "T": ["ACU", "ACC", "ACA", "ACG"],
        "A": ["GCU", "GCC", "GCA", "GCG"],
        "Y": ["UAU", "UAC"],
        "*": ["UAA", "UAG", "UGA"],
        "H": ["CAU", "CAC"],
        "Q": ["CAA", "CAG"],
        "N": ["AAU", "AAC"],
        "K": ["AAA", "AAG"],
        "D": ["GAU", "GAC"],
        "E": ["GAA", "GAG"],
        "C": ["UGU", "UGC"],
        "W": ["UGG"],
        "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "G": ["GGU", "GGC", "GGA", "GGG"],
    }
}


class GeneticCode:
    def __init__(self, name: str = "standard"):
        self._gencode = GENCODE[name]

    def codons(self, amino_acid: str) -> List[str]:
        amino_acid = amino_acid.upper()
        return self._gencode.get(amino_acid, [])


def generate_codon_lprobs(aa_lprobs: Dict[str, float], gencode: GeneticCode):
    from math import log

    codon_lprobs: Dict[str, float] = {}
    for aa, logp in aa_lprobs.items():
        codons = gencode.codons(aa)
        if len(codons) == 0:
            continue
        logp_norm = log(len(codons))
        for codon in codons:
            codon_lprobs[codon] = logp - logp_norm

    return codon_lprobs
