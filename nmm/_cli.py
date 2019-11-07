import click


@click.group(name="nmm", context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    pass


@click.command()
@click.argument("profile", type=click.File("r"))
@click.argument("target", type=click.File("r"))
@click.option("--epsilon", type=float, default=1e-2)
def match(profile, target, epsilon: float):
    """
    Match nucleotide sequences against a HMMER3 Protein profile.
    """
    from nmm import create_frame_profile, read_hmmer
    from fasta_reader import open_fasta

    prof = create_frame_profile(read_hmmer(profile), epsilon=epsilon)

    with open_fasta(target) as fasta:
        for item in fasta:
            r = prof.lr(item.sequence.encode())
            frags = r.fragments

            print(item.defline)
            print(item.sequence)
            print(f"Fragments: {frags}")

            for n, frag in enumerate(frags):
                print(f"Fragment {n}")
                states = []
                matches = []
                for i in frag.items():
                    states.append(i[1].state.name.decode())
                    matches.append(i[0].decode())

                print("\t".join(states))
                print("\t".join(matches))


cli.add_command(match)
