'''Test SNPs specifically in the first and last nucleotides of the
reference sequences.
'''
import common as c
import substitutions as subs


def possible_nt_substitutions(seq, pos):
    org_nt = seq[pos]
    return [nt for nt in c.NUCLEOTIDES if nt != org_nt.upper()]


def possible_mutations(genotype, gene, gene_seq, pos):
    for nt in possible_nt_substitutions(gene_seq, pos):
        orig_codon = subs.codon_at(pos, gene_seq)
        sub_codon = orig_codon[:pos % 3] + nt + orig_codon[(pos % 3) + 1:]
        yield subs.Mutation(
            nt_pos=0,
            aa_pos=0,
            org_nt=gene_seq[pos],
            org_cod=gene_seq[:3],
            sub_nt=nt,
            sub_cod=sub_codon,
            gene=gene,
            genotype=genotype,
        )


class TestSNPsAtBoundaries(c.TestCaseWithReferenceSeqs):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.genotype_seqs = {}
        for seq in cls.reference_seqs:
            genotype = seq.description.split(' ')[-1]
            cls.genotype_seqs[genotype] = seq

    def possible_seqs(self, genotype, gene, pos):
        bio_seq = self.genotype_seqs[genotype]
        seq = subs.gene_seq(bio_seq.seq, genotype, gene)
        for mtn in possible_mutations(genotype, gene, seq, pos):
            yield (mtn, subs.apply_mutation(mtn, seq))

    def check_snp_at(self, genotype, gene, pos):
        for mtn, mtd in self.possible_seqs(genotype, gene, pos):
            profile = 'hcv' + genotype.lower()
            if type(mtd) is not str:
                mtd = str(mtd)
            alignment = c.NucAlignment(mtd, gene, profile)
            self.assertEqual(len(alignment.mutations()), 1)
            report = alignment.mutations()[0]
            mutations = report['mutations']
            if len(mutations) == 0 and alignment.report['FirstAA'] == 2:
                continue
            self.assertEqual(
                len(mutations),
                1,
                "Expected a single mutation",
            )
            mutation = mutations[0]
            self.assertFalse(
                mutation['IsInsertion'] or mutation['IsDeletion'],
                "Unexpected InDel",
            )
            self.assertEqual(mtn.sub_cod, mutation['CodonText'])

    def test_snps_at_beginning(self):
        failures = []
        for genotype in c.GENOTYPES:
            for gene in c.GENES:
                try:
                    self.check_snp_at(genotype, gene, 0)
                except AssertionError as err:
                    failures.append((err, genotype, gene))
        if len(failures) > 1:
            self.fail("Too many snp at beginning failures")
