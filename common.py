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

    @classmethod
    def _nucalign(cls, inputseq, profile='hcv1a', gene="NS3"):
        command = ["./nucamino", "align", profile, gene, "-q", "-f", "json"]
        if type(inputseq) is not bytes:
            inputseq = bytes(inputseq, 'utf8')
        align_proc = subprocess.run(
            command,
            input=inputseq,
            stdout=subprocess.PIPE,
        )
        outp = align_proc.stdout.decode('utf8')
        return json.loads(outp)

    def __init__(self, seq, gene, profile):
        self.nuc_result = self._nucalign(seq, gene=gene, profile=profile)

    def mutations(self):
        for gene, results in self.nuc_result.items():
            mtns = [mtn for r in results for mtn in r['Report']['Mutations']]
            yield {
                'gene': gene,
                'mutations': mtns,
            }
