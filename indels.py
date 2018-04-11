import collections
import random

import Bio.Seq as bioseq

import common

# Use a Tuple beccause Set doesn't support indexing.
_nucleotides = tuple(common.NUCLEOTIDES)


def random_nt():
    return random.choice(_nucleotides)


class BaseIndel(object):
    '''Common properties of insertions and deletions'''

    @property
    def aa_pos(self):
        start, _ = common.GENE_POS[self.genotype][self.gene]
        return (self.nt_pos - start) // 3

    @property
    def delta_len(self):
        raise NotImplementedError()

    @property
    def reference_sequence(self):
        return common.REFSEQS[self.genotype]

    @property
    def gene_pos(self):
        return common.GENE_POS[self.genotype][self.gene]

    @property
    def reference_gene(self):
        start, end = self.gene_pos
        return self.reference_sequence[start:end]

    @property
    def mutated_sequence(self):
        ins = getattr(self, "nt_ins", "")
        skip = getattr(self, "nt_count", 0)
        seq = self.reference_sequence
        pos = self.nt_pos
        return seq[:pos] + ins + seq[pos+skip:]

    @property
    def mutated_gene(self):
        start, end = self.gene_pos
        return self.mutated_sequence[start:end+self.delta_len]


_insertion = collections.namedtuple(
    "Insertion",
    [
        "nt_ins",
        "nt_pos",
        "gene",
        "genotype"
    ],
)


class Insertion(_insertion, BaseIndel):

    def __str__(self):
        return "Ins: {gene} {pos}>{ins} (AA= {aapos})".format(
            gene=self.gene,
            pos=self.nt_pos,
            ins=self.nt_ins,
            aapos=self.aa_pos,
        )

    @classmethod
    def _random(cls, gene, genotype, max_length=1):
        start, end = common.GENE_POS[genotype][gene]
        # Only check whole-codon insertions
        ins_length = 3 * random.randint(1, max_length)
        pos_candidate = start + random.randint(
            1 + 6,
            end - start + ins_length - 6,
        )  # Keep indels at least 2 codons away from the ends.

        def on_codon_boundary(idx):
            return (idx - start) % 3 == 0

        pos = next(idx for idx in range(pos_candidate, -1, -1)
                   if on_codon_boundary(idx))
        ins = "".join(random_nt() for _ in range(ins_length))
        return cls(
            nt_ins=ins,
            nt_pos=pos,
            gene=gene,
            genotype=genotype,
        )

    @property
    def delta_len(self):
        return len(self.nt_ins)

    @property
    def aa_ins(self):
        return bioseq.translate(self.nt_ins)


_deletion = collections.namedtuple(
    "Deletion",
    [
        "nt_pos",
        "gene",
        "genotype",
        "nt_count",
        "orig_nt",
    ],
)


class Deletion(_deletion, BaseIndel):

    # Distance from gene boundaries where random deletions won't
    # generate.
    buffer = 12

    def __str__(self):
        return "Del: {gene} {orig_nt}{pos}{dels} (AA= {aapos})".format(
            gene=self.gene,
            orig_nt=self.orig_nt,
            pos=self.nt_pos,
            dels="-" * self.nt_count,
            aapos=self.aa_pos,
        )

    @property
    def delta_len(self):
        return -1 * self.nt_count

    @classmethod
    def _random(cls, gene, genotype, max_length=1):
        start, end = common.GENE_POS[genotype][gene]
        seq = common.REFSEQS[genotype]
        del_length = 3 * random.randint(1, max_length)
        pos_candidate = start + random.randint(
            1 + cls.buffer,
            end - start - del_length - cls.buffer,
        )

        def on_codon_boundary(idx):
            return (idx - start) % 3 == 0

        pos = next(idx for idx in range(pos_candidate, -1, -1)
                   if on_codon_boundary(idx))
        orig_nt = seq[pos:pos+del_length]
        return cls(
            nt_pos=pos,
            gene=gene,
            genotype=genotype,
            nt_count=del_length,
            orig_nt=orig_nt,
        )

    @property
    def orig_aa(self):
        return bioseq.translate(self.orig_nt)
