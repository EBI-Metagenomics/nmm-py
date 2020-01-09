from typing import Any, Union

import click
from click.utils import LazyFile

from fasta_reader import FASTAWriter
from nmm._gencode import GeneticCode
from nmm._gff import GFFItem, GFFWriter
from nmm._hmmer import SearchResult


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option("--epsilon", type=float, default=1e-2)
@click.option("--output", type=click.File("w"))
@click.option("--ocodon", type=click.File("w"))
@click.option("--oamino", type=click.File("w"))
def search(profile, target, epsilon: float, output, ocodon, oamino):
    """
    Search nucleotide sequences against a HMMER3 Protein profile.
    """
    from nmm import create_frame_profile
    from hmmer_reader import open_hmmer

    record_writer = RecordWriter(epsilon)

    codon_writer: Union[FASTAWriter, None] = None
    if ocodon is not None:
        codon_writer = FASTAWriter(ocodon)

    amino_writer: Union[FASTAWriter, None] = None
    if oamino is not None:
        amino_writer = FASTAWriter(oamino)

    gcode = GeneticCode()

    with open_hmmer(profile) as hmmfile:
        for hmmprof in hmmfile:
            record_writer.accession = hmmprof.metadata["ACC"]
            prof = create_frame_profile(hmmprof, epsilon=epsilon)

            show_header1("Profile")
            print()
            print(hmmprof)
            print()

            show_header1("Targets")

            process_sequence(
                prof, target, record_writer, codon_writer, amino_writer, gcode
            )

            print()

    if output is not None:
        record_writer.dump(output)
        finalize_stream(output)

    if codon_writer is not None:
        finalize_stream(ocodon)

    if amino_writer is not None:
        finalize_stream(oamino)


cli.add_command(search)


class RecordWriter:
    def __init__(self, epsilon: float):
        self._gff = GFFWriter()
        self._accession = "NOTSET"
        self._epsilon = epsilon

    @property
    def accession(self):
        return self._accession

    @accession.setter
    def accession(self, accession: str):
        self._accession = accession

    def add_items(self, result: SearchResult, seqid: str):
        source = f"nmm:{self.accession}"
        for i, frag in enumerate(result.fragments):
            if not frag.homologous:
                continue

            start = result.intervals[i].start
            end = result.intervals[i].end

            att = f"Epsilon={self._epsilon}"
            item = GFFItem(seqid, source, ".", start + 1, end, 0.0, "+", ".", att)
            self._gff.append(item)

    def dump(self, fp):
        self._gff.dump(fp)


class TargetWriter:
    def __init__(self, ocodon=Union[LazyFile, Any], oamino=Union[LazyFile, Any]):
        self._ocodon = ocodon
        self._oamino = oamino

    @property
    def ocodon(self):
        return self._ocodon

    @property
    def oamino(self):
        return self._oamino

    @property
    def has_ocodon(self):
        return self._ocodon is not None

    @property
    def has_oamino(self):
        return self._oamino is not None

    def add_ocodon_target(self, seqid: str, sequence: str):
        if self._ocodon is None:
            return

        self._ocodon.write(">" + seqid + "_codon\n")
        self._ocodon.write(sequence + "\n")

    def add_oamino_target(self, seqid: str, sequence: str):
        if self._oamino is None:
            return

        self._oamino.write(">" + seqid + "_amino\n")
        self._oamino.write(sequence + "\n")


def process_sequence(
    prof,
    target,
    record: RecordWriter,
    codon_writer: Union[FASTAWriter, None],
    amino_writer: Union[FASTAWriter, None],
    gcode,
):
    from fasta_reader import open_fasta

    with open_fasta(target) as fasta:
        for ti, target in enumerate(fasta):
            print()
            show_header2(f"Target {ti}")
            print()

            print(">" + target.defline)
            print(sequence_summary(target.sequence))

            seq = target.sequence.encode().replace(b"T", b"U")
            frame_result = prof.search(seq)
            codon_result = frame_result.decode()
            seqid = f"{target.defline.split()[0]}"

            show_search_result(frame_result)
            record.add_items(frame_result, seqid)

            if codon_writer is not None:
                show_search_result(codon_result)
                ident = seqid + "_codon"
                codon_writer.write_item(ident, codon_result.sequence.decode())
                record.add_items(codon_result, ident)

            if amino_writer is not None:
                amino_result = codon_result.decode(gcode)
                show_search_result(amino_result)
                ident = seqid + "_amino"
                amino_writer.write_item(ident, amino_result.sequence.decode())
                record.add_items(amino_result, ident)


def finalize_stream(stream: Union[LazyFile, Any]):
    if not isinstance(stream, LazyFile):
        return

    if stream.name != "-":
        print(f"Writing to <{stream.name}> file.")

    stream.close_intelligently()


def write_target(file, defline: str, sequence: str):
    file.write(">" + defline + "\n")
    file.write(sequence + "\n")


def show_search_result(result: SearchResult):
    frags = result.fragments
    nhomo = sum(frag.homologous for frag in result.fragments)

    print()
    print(f"Found {nhomo} homologous fragments ({len(frags)} in total).")

    for i, frag in enumerate(frags):
        if not frag.homologous:
            continue
        start = result.intervals[i].start
        end = result.intervals[i].end
        print(f"Homologous fragment={i}; Position=[{start + 1}, {end}]")
        states = []
        matches = []
        for subseq, step in frag.items():
            states.append(step.state.name.decode())
            matches.append(subseq.decode())

        print("\t".join(states))
        print("\t".join(matches))


def show_header1(title: str):
    print(title)
    print("=" * len(title))


def show_header2(title: str):
    print(title)
    print("-" * len(title))


def sequence_summary(sequence: str):
    max_nchars = 79
    if len(sequence) <= max_nchars:
        return sequence

    middle = " ... "

    begin_nchars = (max_nchars - len(middle)) // 2
    end_nchars = begin_nchars + (max_nchars - len(middle)) % 2

    return sequence[:begin_nchars] + middle + sequence[-end_nchars:]
