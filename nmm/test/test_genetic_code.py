from nmm import CanonicalAminoAlphabet, Codon, DNAAlphabet, GeneticCode, RNAAlphabet


def test_genetic_code_dna():
    base_abc = DNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    gcode = GeneticCode(base_abc, amino_abc, "standard")

    assert len(gcode.codons(b"P")) == 4

    assert Codon.create(b"CCT", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in gcode.codons(b"P")

    assert gcode.amino_acid(Codon.create(b"ATG", base_abc)) == b"M"
    assert len(gcode.amino_acids()) == 21
    assert b"R" in gcode.amino_acids()
    assert b"*" in gcode.amino_acids()


def test_genetic_code_rna():
    base_abc = RNAAlphabet()
    amino_abc = CanonicalAminoAlphabet()

    gcode = GeneticCode(base_abc, amino_abc, "standard")

    assert len(gcode.codons(b"P")) == 4

    assert Codon.create(b"CCU", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCC", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCA", base_abc) in gcode.codons(b"P")
    assert Codon.create(b"CCG", base_abc) in gcode.codons(b"P")

    assert gcode.amino_acid(Codon.create(b"AUG", base_abc)) == b"M"
    assert len(gcode.amino_acids()) == 21
    assert b"R" in gcode.amino_acids()
    assert b"*" in gcode.amino_acids()
