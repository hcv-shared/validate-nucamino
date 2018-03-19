import common
import unittest

class TestReferenceSequence(unittest.TestCase):

    reference_filename = 'hcv1a.fasta'

    def test_no_mutations_in_reference(self):
        for gene in ['NS3', 'NS5A', 'NS5B']:
            alignment = common.NucAlignment(self.reference_filename, gene=gene)
            reports = list(alignment.mutations())
            for report in reports:
                for mut in list(report['mutations']):
                    self.assertEqual(len(mut), 0)
