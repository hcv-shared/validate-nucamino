import collections
import functools
import os
import random
import unittest

import Bio.Seq


# ---------------------------------------------------------------------
# Mutations (SNP)
# ---------------------------------------------------------------------
# Substitutions

def make_mut_in_file(infile, outfile='hcv1a_mut.fasta'):
    with open(infile, 'r') as inf:
        inf_header, inf_seq = inf.readlines()
    gene = random.choice(GENES)
    mutation = random_mutation(inf_seq, gene=gene)
    idx = mutation.nt_pos
    outf_header = inf_header.strip() + " with mutation: {}\n".format(mutation)
    outf_seq = apply_mutation(mutation, inf_seq)
    with open(outfile, 'w') as outf:
        outf.writelines([outf_header, outf_seq])
    return mutation


class TestSimpleSubstitutions(unittest.TestCase):

    reference_file = "hcv1a.fasta"
    tmp_filename = "hcv1a_mut.fasta"
    iterations = 1000

    def align_simple_sub(self):
        mutation = make_mut_in_file(
            self.reference_file,
            outfile=self.tmp_filename,
        )
        alignment = NucAlignment(self.tmp_filename, gene=mutation.gene)
        return mutation, alignment

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

    @print_seed_on_assertionerror
    def check_simple_sub(self):
        # Put the RNG in a known state so we can reproduce test
        # failures
        # Known problem seeds:
        # 361de016527e1aa70ad1832874493b9c2290b3b78ebe02e36c451673adb26180 (NS5B S 1 * @ 7602)
        # a7fd60be8525abbc23d82e5c0224b4def33cb13aa52b63ac2259fff553167f0b (NS5A S 1 F @ 6258)

        vmtn, aln = self.align_simple_sub()
        if vmtn.org_aa == vmtn.sub_aa:
            self.check_synonymous_report(aln, vmtn.gene)
            return
        try:
            report = self.check_mutation_report(aln, gene=vmtn.gene)
        except AssertionError as e:
            print(list(aln.mutations()))
            print(vmtn)
            print(vmtn.org_aa)
            print(vmtn.sub_aa)
            raise e
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
        for i in range(self.iterations):
            if i % 20 == 0 and i > 0:
                print(i)
            self.check_simple_sub()


# ---------------------------------------------------------------------
# Generate and Work with Indels

_insertion = collections.namedtuple(
    "Insertion",
    [
        "nt_ins",
        "nt_pos",
        "gene",
    ],
)


class Insertion(_insertion):

    @property
    def aa_pos(self):
        start, _ = GENE_POS['1a'][self.gene]
        return (self.nt_pos - start) // 3

    def __str__(self):
        return "{gene} {pos}>{ins} (NA= {napos})".format(
            gene=self.gene,
            pos=self.nt_pos,
            ins=self.nt_ins,
            napos=self.aa_pos,
        )


def random_nt():
    return random.choice("GATC")


def apply_insertion(seq, insertion):
    pos = insertion.nt_pos
    ins = insertion.nt_ins
    return seq[:pos] + ins + seq[pos:]


def apply_deletion(seq, pos, count):
    return seq[:pos] + seq[pos+count:]


def random_insertion(seq, max_length=1, gene=None, genotype='1a'):
    start, end = GENE_POS[genotype][gene]
    ins_length = random.randint(1, max_length)
    pos = start + random.randint(1, (end - start - ins_length))
    ins = "".join(random_nt() for _ in range(ins_length))
    return Insertion(
        nt_ins=ins,
        nt_pos=pos,
        gene=gene,
    )


def random_deletion(seq, max_length=1, gene=None):
    count = random.randint(1, max_length)
    pos = random.choice(range(len(seq) - count))
    return (pos, count)


class TestIndelGeneration(unittest.TestCase):

    def test_apply_insertion(self):
        base_seq = "abc"
        cases = [
            (0, 'x', 'xabc'),
            (1, 'x', 'axbc'),
            (2, 'x', 'abxc'),
            (3, 'x', 'abcx'),
        ]
        for pos, nt_ins, expected in cases:
            insertion = Insertion(nt_ins=nt_ins, nt_pos=pos, gene=None)
            inserted = apply_insertion(base_seq, insertion)
            self.assertEqual(inserted, expected)

    def test_apply_deletion(self):
        base_seq = "abc"
        cases = [
            (0, 1, "bc"),
            (0, 2, "c"),
            (0, 3, ""),
            (1, 1, "ac"),
            (1, 2, "a"),
            (2, 1, "ab"),
        ]
        for pos, count, expected in cases:
            deleted = apply_deletion(base_seq, pos, count)
            self.assertEqual(deleted, expected)

    def check_insertion_generation(self):
        max_len = random.randint(1, 10)
        gene = random.choice(GENES)
        start, end = GENE_POS['1a'][gene]
        insertion = random_insertion(
            "",
            max_length=max_len,
            gene=gene,
        )
        self.assertLessEqual(start, insertion.nt_pos)
        self.assertLessEqual(insertion.nt_pos, end)
        self.assertLessEqual(len(insertion.nt_ins), max_len)
        self.assertTrue(set(insertion.nt_ins).issubset(NUCLEOTIDES))

    @print_seed_on_assertionerror
    def test_insertion_generation(self):
        for i in range(10000):
            self.check_insertion_generation()

    def check_deletion_generation(self):
        max_count = random.rand_int(1, 10)
        gene = random.choice(GENES)
        start, end = GENE_POS['1a'][gene]
        rand_pos, rand_count = random_deletion(
            "",
            max_length=max_length,
            gene=gene,
        )
        self.assertLessEqual(start, rand_pos)
        self.assertLessEqual(rand_pos, end)
        self.assertLessEqual(rand_pos + rand_count, end)
        self.assertLessEqual(rand_count, max_count)

    @print_seed_on_assertionerror
    def test_deletion_generation(self):
        for i in range(10000):
            self.check_insertion_generation()


# ---------------------------------------------------------------------
# Test Alinger on Indels

class TestInsertions(unittest.TestCase):

    tmp_filename = "hcv1a_mut.fasta"
    iterations = 100000

    def setUp(self):
        with open("hcv1a.fasta") as inf:
            _, self.hcv1a_seq = inf.readlines()

    def create_file_with_insertion(self):
        gene = random.choice(GENES)
        ins = random_insertion(
            self.hcv1a_seq,
            max_length=5,
            gene=gene,
        )
        applied = apply_insertion(self.hcv1a_seq, ins)
        with open(self.tmp_filename, "w") as outf:
            outf.writelines([
                "> HCV1A with insertion {}\n".format(ins),
                applied,
            ])
        return ins

    def check_make_insertion_file(self):
        insertion = self.create_file_with_insertion()
        pos = insertion.nt_pos
        ins = insertion.nt_ins
        with open("hcv1a.fasta") as infile:
            _, orig = infile.readlines()
        with open(self.tmp_filename) as infile:
            _, seq = infile.readlines()
        self.assertEqual(
            len(orig) + len(ins),
            len(seq),
        )

        self.assertEqual(ins, seq[pos:pos+len(ins)])

    @print_seed_on_assertionerror
    def test_make_insertion_file(self):
        # for _ in range(self.iterations):
        self.check_make_insertion_file()
