import click


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option("--epsilon", type=float, default=1e-2)
@click.option(
    # "--output", type=click.Path(file_okay=True, dir_okay=False, writable=True)
    "--output",
    type=click.File("w"),
)
def match(profile, target, epsilon: float, output):
    """
    Match nucleotide sequences against a HMMER3 Protein profile.
    """
    from nmm import create_frame_profile, read_hmmer
    from nmm._gff import GFFWriter, Item as GFFItem
    from fasta_reader import open_fasta

    hmmer_reader = read_hmmer(profile)
    prof_acc = hmmer_reader.metadata["ACC"]
    prof = create_frame_profile(hmmer_reader, epsilon=epsilon)

    print("Profile")
    print("=======")
    print()
    print(hmmer_reader)

    print()
    print("Targets")
    print("=======")

    gff = GFFWriter()
    with open_fasta(target) as fasta:
        for ti, target in enumerate(fasta):
            print()
            section = f"Target {ti}"
            print(section)
            print("-" * len(section))
            print()

            print(">" + target.defline)
            print(target.sequence)

            r = prof.lr(target.sequence.encode())
            frags = r.fragments

            hfrags = [frag for frag in r.fragments if frag.homologous]

            print()
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

                gff.append(
                    GFFItem(
                        seqid=f"{target.defline}",
                        source=f"nmm:{prof_acc}",
                        type=".",
                        start=start,
                        end=end,
                        score=0.0,
                        strand="+",
                        phase=".",
                        attributes=".",
                    )
                )
    gff.dump(output)


cli.add_command(match)
