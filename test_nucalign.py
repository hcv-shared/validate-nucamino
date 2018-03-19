import collections
import functools
import os
import random
import unittest

import Bio.Seq

import common


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
