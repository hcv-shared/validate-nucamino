import random

import Bio.Seq as bioseq

import common
import indels


class BaseTestIndel(common.TestCaseWithReferenceSeqs):

    iterations = 10

    def check_aligning_indels(self, label):
        failures = []
        for genotype in common.GENOTYPES:
            print("Testing {}".format(genotype), end="", flush=True)
            for i in range(self.iterations):
                if i % 20 == 0 and i > 0:
                    print(".", end="", flush=True)
                mtn, aln = self.align(genotype)
                try:
                    self.check_alignment(mtn, aln)
                except AssertionError:
                    failures.append((genotype, mtn, aln))
            print(flush=True)
        self.check_genotype_failures(failures, label=label)

    def get_unique_aligned_mutation(self, mutation, alignment):
        mtn_reports = alignment.mutations()
        self.assertEqual(len(mtn_reports), 1, "Expected one report")
        mtn_report = mtn_reports[0]
        self.assertEqual(mutation.gene, mtn_report["gene"])
        mtns = mtn_report["mutations"]
        self.assertEqual(len(mtns), 1, "Expected one mutation")
        mtn = mtns[0]
        return mtn


class TestInsertions(BaseTestIndel):

    maxDiff = 16_000   # Display sequence differences

    def align(self, genotype):
        gene = random.choice(common.GENES)
        insertion = indels.Insertion._random(
            gene,
            genotype,
            max_length=2,
        )
        mtd = insertion.mutated_sequence
        profile = "hcv" + genotype
        alignment = common.NucAlignment(mtd, gene=gene, profile=profile)
        return insertion, alignment

    def check_alignment(self, insertion, alignment):
        mtn = self.get_unique_aligned_mutation(insertion, alignment)
        try:
            self.insertion_is_exact_match(insertion, mtn)
        except AssertionError:
            self.insertion_is_synonymous_match(insertion, mtn)

    def insertion_is_exact_match(self, insertion, mtn):
        self.assertEqual(
            len(mtn["InsertedCodonsText"]),
            len(insertion.nt_ins),
            "Virtual and aligned insertion lengths don't match",
        )
        # To determine if the mutation nucleotide position is the same
        # as the aligned mutation's position, increment it by one to
        # account for 0/1 indexing, and decrease it by 3 because
        # nucamino reports the insertion as beginning at the reference
        # sequence's previous codon.
        self.assertEqual(
            insertion.nt_pos + 1 - 3,
            mtn['NAPosition'],
            "Position mismatch",
        )
        self.assertEqual(
            insertion.nt_ins,
            mtn['InsertedCodonsText'],
            "Generated insertion ≠ Aligned mutation",
        )

    def insertion_is_synonymous_match(self, insertion, mtn):
        nt_pos = mtn["NAPosition"] + 3 - 1
        nt_ins = mtn["InsertedCodonsText"]
        recovered_insertion = indels.Insertion(
            nt_ins=nt_ins,
            nt_pos=nt_pos,
            gene=insertion.gene,
            genotype=insertion.genotype,
        )
        self.assertGreater(len(nt_ins), 0)
        mutated = insertion.mutated_gene
        aligned = recovered_insertion.mutated_gene
        self.assertEqual(
            bioseq.translate(mutated),
            bioseq.translate(aligned),
        )


    def test_aligning_insertions(self):
        self.check_aligning_indels(label="ins")


class TestDeletions(BaseTestIndel):

    def align(self, genotype):
        gene = random.choice(common.GENES)
        deletion = indels.Deletion._random(gene, genotype, max_length=2)
        mtd = deletion.mutated_sequence
        profile = "hcv" + genotype
        alignment = common.NucAlignment(mtd, gene=gene, profile=profile)
        return deletion, alignment

    def check_alignment(self, deletion, alignment):
        mtn = self.get_unique_aligned_mutation(deletion, alignment)
        self.assertTrue(mtn["IsDeletion"], "Expected a deletion")
        try:
            self.deletion_is_exact_match(deletion, mtn)
        except AssertionError:
            self.deletion_is_synonymous_match(deletion, mtn)

    def deletion_is_exact_match(self, deletion, mtn):
        self.assertEqual(
            deletion.nt_pos + 1,  # ∵ nucamino is 1-indexed
            mtn["NAPosition"],
            "Position mismatch",
        )
        self.assertEqual(
            deletion.nt_count,
            mtn["Control"].count('-'),
            "Length mismatch",
        )

    def deletion_is_synonymous_match(self, deletion, mtn):
        nt_pos = mtn["NAPosition"] - 1
        nt_count = mtn["Control"].count('-')
        recovered_deletion = indels.Deletion(
            nt_pos=nt_pos,
            gene=deletion.gene,
            genotype=deletion.genotype,
            nt_count=nt_count,
            orig_nt=None,
        )
        self.assertGreater(nt_count, 0)
        mutated = deletion.mutated_gene
        aligned = recovered_deletion.mutated_gene
        # Translate to amino-acids to account for synonymous mutations.
        self.assertEqual(
            bioseq.translate(mutated),
            bioseq.translate(aligned),
        )

    def test_aligning_deletions(self):
        self.check_aligning_indels(label="dels")
