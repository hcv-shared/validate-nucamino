import collections
import random
import unittest

import common

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
        start, _ = common.GENE_POS['1a'][self.gene]
        return (self.nt_pos - start) // 3

    def __str__(self):
        return "{gene} {pos}>{ins} (NA= {napos})".format(
            gene=self.gene,
            pos=self.nt_pos,
            ins=self.nt_ins,
            napos=self.aa_pos,
        )


_nucleotides = tuple(common.NUCLEOTIDES)


def random_nt():
    return random.choice(_nucleotides)


def apply_insertion(seq, insertion):
    pos = insertion.nt_pos
    ins = insertion.nt_ins
    return seq[:pos] + ins + seq[pos:]


def apply_deletion(seq, pos, count):
    return seq[:pos] + seq[pos+count:]


def random_insertion(seq, max_length=1, gene=None, genotype='1a'):
    start, end = common.GENE_POS[genotype][gene]
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
        gene = random.choice(common.GENES)
        start, end = common.GENE_POS['1a'][gene]
        insertion = random_insertion(
            "",
            max_length=max_len,
            gene=gene,
        )
        self.assertLessEqual(start, insertion.nt_pos)
        self.assertLessEqual(insertion.nt_pos, end)
        self.assertLessEqual(len(insertion.nt_ins), max_len)
        self.assertTrue(set(insertion.nt_ins).issubset(common.NUCLEOTIDES))

    @common.print_seed_on_assertionerror
    def test_insertion_generation(self):
        for i in range(10000):
            self.check_insertion_generation()

    def check_deletion_generation(self):
        max_count = random.rand_int(1, 10)
        gene = random.choice(common.GENES)
        start, end = common.GENE_POS['1a'][gene]
        rand_pos, rand_count = random_deletion(
            "",
            max_length=max_length,
            gene=gene,
        )
        self.assertLessEqual(start, rand_pos)
        self.assertLessEqual(rand_pos, end)
        self.assertLessEqual(rand_pos + rand_count, end)
        self.assertLessEqual(rand_count, max_count)

    @common.print_seed_on_assertionerror
    def test_deletion_generation(self):
        for i in range(10000):
            self.check_insertion_generation()
