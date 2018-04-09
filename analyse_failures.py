import collections
import pickle

import Bio.SeqIO as seqio
import Bio.Seq as bioseq

import common
import indels
import test_indels


def load_failures(fname="failures.pickle"):
    with open(fname, 'rb') as inf:
        failures = pickle.load(inf)
    return failures


def load_refseqs():
    with open("hcv-refs.fasta") as inf:
        seqs = list(seqio.parse(inf, 'fasta'))
    result = {}
    for seq in seqs:
        key = seq.description.split(' ')[-1]
        seq = str(seq.seq)
        result[key] = seq
    return result


class FailData(object):

    def __init__(self, failure):
        self.genotype, self.mutation, self.alignment = failure

    @property
    def result(self):
        nuc_result = self.alignment.nuc_result
        assert len(nuc_result) == 1
        nuc_result = next(iter(nuc_result.values()))
        assert len(nuc_result) == 1
        result = next(iter(nuc_result))
        return result

    @property
    def report(self):
        return self.result['Report']

    @property
    def aln_muts(self):
        return self.report["Mutations"]

    @property
    def aln_pos(self):
        muts = self.aln_muts
        if len(muts) != 1:
            return None
        return muts[0]["NAPosition"]

    @property
    def mut_pos(self):
        self.mutation.nt_pos + 1  # ∵ 1 indexing

    @property
    def distance(self):
        return self.mut_pos - self.aln_pos

    @property
    def muts_and_frameshifts(self):
        muts = len(self.report['Mutations'])
        fs = len(self.report['FrameShifts'])
        return muts, fs

    @property
    def blank(self):
        return self.muts_and_frameshifts == (0, 0)

    @property
    def refseq(self):
        return refseqs[self.mutation.genotype]

    @property
    def gene_pos(self):
        genotype = self.mutation.genotype
        gene = self.mutation.gene
        start, end = common.GENE_POS[genotype][gene]
        return (start, end)

    @property
    def fasta_text(self):
        start, end = self.gene_pos
        ref = self.refseq
        mutated = self.mutation.mutated_sequence
        aligned = self.report['NucleicAcidsLine']
        parts = [
            "> Reference ({})\n".format(self.mutation.genotype),
            ref[start:end+1],
            "> Mutated ({})\n".format(self.mutation),
            mutated[start:(end+1)+self.mutation.delta_len],
            "> Aligned",
            aligned,
        ]
        return "\n".join(parts)

    def save_fasta_text(self, fname="failure.fasta"):
        with open(fname, 'w') as outf:
            outf.write(self.fasta_text)


failures = load_failures("failures-dels.pickle")
refseqs = load_refseqs()

data = [FailData(f) for f in failures]
d = data[0]
assert len(collections.Counter(d.result['Error'] for d in data)) == 1
assert len(collections.Counter(d.result['Err'] for d in data)) == 1
assert sum(len(d.report['FrameShifts']) for d in data) == 0
