import random
import unittest

import Bio.SeqIO as seqio

import common
import substitutions


def random_mut_in_seq(seq, genotype):
    gene = random.choice(common.GENES)
    mutation = substitutions.random_mutation(seq, gene, genotype)
    mutated = substitutions.apply_mutation(mutation, seq)
    return mutation, mutated


def align_simple_sub(seq, genotype):
    mtn, mtd = random_mut_in_seq(seq, genotype)
    profile = "hcv" + mtn.genotype
    alignment = common.NucAlignment(mtd, gene=mtn.gene, profile=profile)
    return mtn, alignment


class TestSimpleSubstitutions(unittest.TestCase):

    reference_file = "hcv-refs.fasta"
    iterations = 400

    @classmethod
    def setUpClass(cls):
        with open(cls.reference_file) as inf:
            seqs = list(seqio.parse(inf, 'fasta'))
        cls.reference_seqs = seqs

    def check_vmtn_mtn(self, vmtn, mtn):
        self.assertEqual(vmtn.aa_pos, mtn['Position'])
        self.assertEqual(vmtn.org_aa, mtn['ReferenceText'])
        self.assertEqual(vmtn.sub_aa, mtn['AminoAcidText'])
        self.assertIn(
            vmtn.nt_pos - mtn['NAPosition'],
            {-1, 0, 1},
        )

    def check_mutation_report(self, alignment, gene=None):
        reports = list(
            r for r in alignment.mutations()
            if len(r['mutations']) > 0
        )
        self.assertEqual(len(reports), 1, "Expected a single mutation report")
        report = next(iter(reports))
        self.assertEqual(
            len(report['mutations']),
            1,
            "Expected a single mutation",
        )
        self.assertEqual(gene, report['gene'])
        return report

    def check_synonymous_report(self, alignment, gene=None):
        reports = list(alignment.mutations())
        for rep in reports:
            msg = "Unexpected mutation in {}".format(rep)
            self.assertFalse(rep['mutations'], msg)

    @common.print_seed_on_assertionerror
    def check_simple_sub(self, seq, genotype):
        vmtn, aln = align_simple_sub(seq, genotype)
        if vmtn.org_aa == vmtn.sub_aa:
            self.check_synonymous_report(aln, vmtn.gene)
            return
        report = self.check_mutation_report(aln, gene=vmtn.gene)
        mtn = next(iter(report['mutations']))
        try:
            self.check_vmtn_mtn(vmtn, mtn)
        except AssertionError as e:
            print(e)
            print(mtn)
            print(aln)
            raise e

    def test_simple_subs(self):
        print()
        for bioseq in self.reference_seqs:
            print("Testing {}".format(bioseq.description))
            seq = str(bioseq.seq)
            genotype = bioseq.description.split()[1]
            for i in range(self.iterations):
                if i % 20 == 0 and i > 0:
                    print(i)
                    self.check_simple_sub(seq, genotype)
