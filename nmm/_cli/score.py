from contextlib import contextmanager

import click


@click.command()
@click.argument("profile", type=click.Path(exists=True))
@click.argument("target", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
def score(profile, target, gff):
    """
    Score sequences against a database of HMMER3 profiles.
    """

    hmmsearch = HMMSearch()
    scores = hmmsearch.compute_scores(profile, target)
    hmmsearch.update_gff_file(gff, scores)


class HMMSearch:
    def __init__(self):
        import distutils.spawn

        program = "hmmsearch"
        prog_path = distutils.spawn.find_executable(program)
        if prog_path is None:
            raise click.ClickException(f"Could not find the `{program}` program.")

        self._prog_path = prog_path

    def compute_scores(self, profile, target):
        import os
        import subprocess
        from nmm import tblout_reader

        profile = os.path.abspath(profile)
        target = os.path.abspath(target)

        with tmp_cwd():
            cmd = [self._prog_path, "--tblout", "tblout", profile, target]
            subprocess.check_output(cmd)
            scores = {}
            with open("tblout", "r") as file:
                for row in tblout_reader(file):
                    scores[row.target_name] = row.full_sequence.e_value

        return scores

    def update_gff_file(self, filepath, scores):
        import fileinput
        import re

        for row in fileinput.input(filepath, inplace=True, backup=".bak"):

            if row.startswith("#"):
                print(row, end="")
                continue

            match = re.search(r"ID=([^;]+);", row)
            if match is None:
                print(row, end="")
                continue

            item_id = match.group(1)
            if item_id not in scores:
                print(row, end="")
                continue

            row = row.rstrip()
            print(row + ";E-value={}".format(scores[item_id]))


@contextmanager
def tmp_cwd():
    import os
    from tempfile import TemporaryDirectory

    oldpwd = os.getcwd()
    with TemporaryDirectory() as tmpdir:

        os.chdir(tmpdir)
        try:
            yield
        finally:
            os.chdir(oldpwd)
