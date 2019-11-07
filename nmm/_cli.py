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
    from nmm import create_frame_profile, read_hmmer, FASTAReader

    prof = create_frame_profile(read_hmmer(profile), epsilon=epsilon)

    with FASTAReader(target) as fasta:
        for item in fasta.items():
            defline = item[0]
            seq = item[1]

            r = prof.lr(seq.encode())
            frags = r.fragments

            print(defline)
            print(seq)
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
