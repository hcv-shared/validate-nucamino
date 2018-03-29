import collections
import functools
import json
import random
import secrets
import subprocess
import unittest

import Bio.SeqIO as seqio


GENOTYPES = ['1a', '1b', '2', '3', '4', '5', '6']
GENES = ["NS3", "NS5A", "NS5B"]
GENE_POS = {
    "1a": {"NS3": (3419, 5312), "NS5A": (6257, 7601), "NS5B": (7601, 9377)},
    "1b": {"NS3": (3419, 5312), "NS5A": (6257, 7598), "NS5B": (7598, 9374)},
    "2": {"NS3": (3430, 5323), "NS5A": (6268, 7666), "NS5B": (7666, 9442)},
    "3": {"NS3": (3435, 5328), "NS5A": (6273, 7629), "NS5B": (7629, 9402)},
    "4": {"NS3": (3418, 5311), "NS5A": (6256, 7591), "NS5B": (7591, 9364)},
    "5": {"NS3": (3327, 5220), "NS5A": (6165, 7515), "NS5B": (7515, 9291)},
    "6": {"NS3": (3373, 5266), "NS5A": (6211, 7564), "NS5B": (7564, 9340)},
}
NUCLEOTIDES = {'G', 'C', 'T', 'A'}


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
        return json.loads(outp), align_proc

    def __init__(self, seq, gene, profile):
        self.nuc_result, self.nuc_proc = self._nucalign(
            seq,
            gene=gene,
            profile=profile,
        )

    def mutations(self):
        for gene, results in self.nuc_result.items():
            mtns = [mtn for r in results for mtn in r['Report']['Mutations']]
            yield {
                'gene': gene,
                'mutations': mtns,
            }


class TestCaseWithReferenceSeqs(unittest.TestCase):

    reference_file = "hcv-refs.fasta"

    @classmethod
    def setUpClass(cls):
        with open(cls.reference_file) as inf:
            seqs = list(seqio.parse(inf, 'fasta'))
        cls.reference_seqs = seqs
        super().setUpClass()

    def check_genotype_failures(self, failures):
        trial_count = len(self.reference_seqs) * self.iterations
        if len(failures) / trial_count > 0.001:
            msg = "More than 0.1% of trials failed ({} / {}".format(
                len(failures),
                trial_count,
            )
            self.fail()
        by_gt = collections.Counter(gt for gt, _ in failures)
        if len(failures) > 0:
            print("Failure rates:")
            for gt, failures in by_gt.items():
                print(f"{gt}  {failures} / {self.iterations}")
            worst_gt, worst_gt_failures = by_gt.most_common()[0]
            if worst_gt_failures / self.iterations > 0.01:
                msg = "{} failed {} / {} of its trials".format(
                    worst_gt,
                    worst_gt_failures,
                    self.iterations,
                )
                self.fail(msg)
