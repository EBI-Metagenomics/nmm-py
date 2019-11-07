import click


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option("--epsilon", type=float, default=1e-2)
@click.option(
    "--output", type=click.Path(file_okay=True, dir_okay=False, writable=True)
)
def match(profile, target, epsilon: float, output):
    """
    Match nucleotide sequences against a HMMER3 Protein profile.
    """
    from nmm import create_frame_profile, read_hmmer
    from nmm._gff import GFFWriter, Item as GFFItem
    from fasta_reader import open_fasta

    prof = create_frame_profile(read_hmmer(profile), epsilon=epsilon)

    gff = GFFWriter()
    with open_fasta(target) as fasta:
        for ti, target in enumerate(fasta):
            print(f"Target: {ti}")
            print(">" + target.defline)
            print(target.sequence)

            r = prof.lr(target.sequence.encode())
            frags = r.fragments

            hfrags = [frag for frag in r.fragments if frag.homologous]

            print(f"Found {len(hfrags)} homologous fragments ({len(frags)} in total).")

            for fi, frag in enumerate(hfrags):
                start = frag.interval.start
                end = frag.interval.end
                print(f"Homologous fragment={fi}; Position=[{start}, {end}]")
                states = []
                matches = []
                for i in frag.items():
                    states.append(i[1].state.name.decode())
                    matches.append(i[0].decode())

                print("\t".join(states))
                print("\t".join(matches))

                gff.append(GFFItem(seqid=f"{target.defline}"))

            print()


cli.add_command(match)
