from typing import Dict, List, Set, Union

from ._alphabet import AminoAlphabet, DNAAlphabet, RNAAlphabet
from ._codon import Codon

__all__ = ["GeneticCode"]


class GeneticCode:
    """
    Genetic code.

    Parameters
    ----------
    base_abc
        Base alphabet.
    amino_abc
        Amino acid alphabet.
    name
        NCBI `translation table name`_. Defaults to `"Standard"`.

    .. _translation table name: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """

    def __init__(
        self,
        base_abc: Union[DNAAlphabet, RNAAlphabet],
        amino_abc: AminoAlphabet,
        name: str = "Standard",
    ):

        self._base_alphabet = base_abc
        self._amino_alphabet = amino_abc

        table = get_translation_table(name)

        self._gencode: Dict[bytes, List[Codon]] = {aa: [] for aa in table}

        for aa, triplets in table.items():
            gcode = self._gencode[aa]
            for triplet in triplets:
                if isinstance(base_abc, RNAAlphabet):
                    triplet = triplet.replace(b"T", b"U")
                gcode.append(Codon.create(triplet, base_abc))

        self._amino_acid: Dict[Codon, bytes] = {}
        for aa, codons in self._gencode.items():
            for codon in codons:
                self._amino_acid[codon] = aa

    def codons(self, amino_acid: bytes) -> List[Codon]:
        amino_acid = amino_acid.upper()
        return self._gencode.get(amino_acid, [])

    def codons_prob(self, amino_acid: bytes) -> Dict[Codon, float]:
        codons = self.codons(amino_acid)
        n = len(codons)
        if n == 0:
            return {}
        return {codon: 1 / n for codon in codons}

    def amino_acid(self, codon: Codon) -> bytes:
        return self._amino_acid[codon]

    def amino_acids(self) -> Set[bytes]:
        return set(self._gencode.keys())

    @property
    def base_alphabet(self) -> Union[DNAAlphabet, RNAAlphabet]:
        return self._base_alphabet

    @property
    def amino_alphabet(self) -> AminoAlphabet:
        return self._amino_alphabet


def get_translation_table(name: str = "Standard"):
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PendingDeprecationWarning)
        from Bio.Data import CodonTable

    if name not in CodonTable.unambiguous_dna_by_name:
        names = str(list(CodonTable.unambiguous_dna_by_name.keys()))
        msg = f"Unknown translation table {name}. Possible names are: {names}."
        raise ValueError(msg)

    table = CodonTable.unambiguous_dna_by_name[name]

    btable: Dict[bytes, List[bytes]] = {}

    for codon, aa in table.forward_table.items():
        baa = aa.encode()
        if baa not in btable:
            btable[baa] = []
        btable[baa].append(codon.encode())

    return btable
