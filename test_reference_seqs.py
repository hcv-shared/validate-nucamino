import common
import unittest

import Bio.SeqIO as seqio

class TestReferenceSequence(unittest.TestCase):

    references_filename = 'hcv-refs.fasta'
    genes = ('NS3', 'NS5A', 'NS5B')

    def test_no_mutations_in_reference(self):
        gene = ','.join(self.genes)
        with open(self.references_filename) as inf:
            refseqs = seqio.parse(inf, 'fasta')
            for seq in refseqs:
                genotype = seq.description.split()[1]
                profile = 'hcv' + genotype
                alignment = common.NucAlignment(
                    str(seq.seq),
                    profile=profile,
                    gene=gene,
                )
                reports = list(alignment.mutations())
                for report in reports:
                    muts = list(report['mutations'])
                    msg = "Found mutations in the reference {}".format(genotype)
                    self.assertEqual(len(muts), 0, msg)
