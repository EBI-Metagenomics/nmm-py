import click
from click.utils import LazyFile
from typing import Union, Any

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
    from nmm import create_frame_profile, read_hmmer
    from fasta_reader import open_fasta

    hmmer_reader = read_hmmer(profile)
    acc = hmmer_reader.metadata["ACC"]
    prof = create_frame_profile(hmmer_reader, epsilon=epsilon)

    show_header1("Profile")
    print()
    print(hmmer_reader)
    print()

    show_header1("Targets")

    gff = GFFWriter()
    gcode = GeneticCode()
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
            create_gffitems(gff, frame_result, seqid, acc, epsilon)

            if ocodon is not None:
                show_search_result(codon_result)
                write_target(ocodon, seqid + "_codon", codon_result.sequence.decode())
                create_gffitems(gff, codon_result, seqid + "_codon", acc, epsilon)

            if oamino is not None:
                amino_result = codon_result.decode(gcode)
                show_search_result(amino_result)
                write_target(oamino, seqid + "_amino", amino_result.sequence.decode())
                create_gffitems(gff, amino_result, seqid + "_amino", acc, epsilon)

    print()
    if output is not None:
        gff.dump(output)
        finalize_stream(output)

    if ocodon is not None:
        finalize_stream(ocodon)

    if oamino is not None:
        finalize_stream(oamino)


cli.add_command(search)


def finalize_stream(stream: Union[LazyFile, Any]):
    if not isinstance(stream, LazyFile):
        return

    if stream.name != "-":
        print(f"Writing to <{stream.name}> file.")

    stream.close_intelligently()


def write_target(file, defline: str, sequence: str):
    file.write(">" + defline + "\n")
    file.write(sequence + "\n")


def create_gffitems(
    gff: GFFWriter, result: SearchResult, seqid: str, accession: str, epsilon: float
):
    source = f"nmm:{accession}"
    for i, frag in enumerate(result.fragments):
        if not frag.homologous:
            continue

        start = result.intervals[i].start
        end = result.intervals[i].end

        att = f"Epsilon={epsilon}"
        item = GFFItem(seqid, source, ".", start + 1, end, 0.0, "+", ".", att)
        gff.append(item)


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
