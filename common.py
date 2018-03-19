import functools
import json
import secrets
import subprocess


def print_seed_on_assertionerror(f):
    "Decorate a test to set the rng seed (and print it if the test fails)"

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        seed = secrets.token_hex()
        random.seed(seed)
        try:
            f(*args, **kwargs)
        except AssertionError as e:
            print("seed = {}".format(seed))
            raise e

    return wrapped


class NucAlignment(object):
    "Perform an alignment on a given sequence and store the results"

    @staticmethod
    def _nucalign(infilename, profile='hcv1a', gene="NS3", output_format="json"):
        command = [
            './nucamino', 'align', profile, gene,
            '-q',
            '--output-format', output_format,
            '-i', infilename,
        ]
        with subprocess.Popen(command, stdout=subprocess.PIPE) as nuc_proc:
            outp = nuc_proc.stdout.read()
            return json.loads(outp)

    def __init__(self, infilename, gene):
        self.nuc_result = self._nucalign(infilename, gene=gene)

    def mutations(self):
        for gene, results in self.nuc_result.items():
            mtns = [mtn for r in results for mtn in r['Report']['Mutations']]
            yield {
                'gene': gene,
                'mutations': mtns,
            }
