import collections
import functools
import json
import os
import random
import secrets
import subprocess
import unittest

import Bio.Seq


class NucAlignment(object):

    @staticmethod
    def nucalign(infilename, cmd='hcv1a', gene="NS3", output_format="json"):
        command = [
            './nucamino', cmd, '-q',
            '-g', gene,
            '--output-format', output_format,
            '-i', infilename,
        ]
        with subprocess.Popen(command, stdout=subprocess.PIPE) as nuc_proc:
            outp = nuc_proc.stdout.read()
            return json.loads(outp)

    def __init__(self, infilename, gene):
        self.nuc_result = self.nucalign(infilename, gene=gene)

    def mutations(self):
        for gene, results in self.nuc_result.items():
            mtns = [mtn for r in results for mtn in r['Report']['Mutations']]
            yield {
                'gene': gene,
                'mutations': mtns,
            }


class TestReferenceSequence(unittest.TestCase):

    reference_filename = 'hcv1a.fasta'

    def test_no_mutations_in_reference(self):
        for gene in ['NS3', 'NS5A', 'NS5B']:
            alignment = NucAlignment(self.reference_filename, gene=gene)
            reports = list(alignment.mutations())
            for report in reports:
                for mut in list(report['mutations']):
                    self.assertEqual(len(mut), 0)


# ---------------------------------------------------------------------

GENE_POS = {
    'NS3': (3419, 5312),
    'NS5A': (6257, 7601),
    'NS5B': (7601, 9374),
}

NUCLEOTIDES = {'G', 'C', 'T', 'A'}


_Mutation = collections.namedtuple(
    "_Mutation",
    ["nt_pos", "aa_pos", "sub_nt", "org_nt", "sub_cod", "org_cod", "gene"],
)

class Mutation(_Mutation):

    @property
    def org_aa(self):
        return Bio.Seq.translate(self.org_cod)

    @property
    def sub_aa(self):
        return Bio.Seq.translate(self.sub_cod)

    def __str__(self):
        sub = "{} {} {} {} @ {}".format(
            self.gene,
            self.org_aa,
            self.aa_pos,
            self.sub_aa,
            self.nt_pos

        )
        return "<Mutation {}>".format(sub)


def gene_seq(seq, gene=None):
    start, end = GENE_POS[gene]
    return seq[start:end]


def random_substitution(orig):
    return random.choice(list(NUCLEOTIDES - {orig}))


def codon_idx(idx):
    return idx - (idx % 3)


def codon_at(idx, seq):
    cod_idx = codon_idx(idx)
    return seq[cod_idx:cod_idx + 3]


def random_mutation(seq, gene='NS3'):
    start, end = GENE_POS[gene]
    gseq = gene_seq(seq, gene=gene)
    idx = random.choice(range(len(gseq)))
    orig_nt = gseq[idx]
    orig_codon = codon_at(idx, gseq)
    sub = random_substitution(orig_nt)
    sub_codon = orig_codon[:idx%3] + sub + orig_codon[(idx%3)+1:]
    return Mutation(
        nt_pos=idx+start,
        aa_pos=idx//3 + 1,
        org_nt=orig_nt,
        org_cod=orig_codon,
        sub_nt=sub,
        sub_cod=sub_codon,
        gene=gene,
    )


def apply_mutation(mut, sequence):
    idx = mut.nt_pos
    return (
        sequence[:idx]
        + mut.sub_nt
        + sequence[idx+1:]
    )


class TestTranslation(unittest.TestCase):

    @staticmethod
    def translate_hcv1a(start, end):
        with open('hcv1a.fasta') as infile:
            _, hcv1a_nt_seq = infile.readlines()
        nt_seq = hcv1a_nt_seq[start:end]
        return Bio.Seq.translate(nt_seq)

    NS3_AA = (
        "APITAYAQQTRGLLGCIITSLTGRDKNQVEGEVQIVSTATQTFLATCINGVCWTVYHGAGTRTIAS"
        "PKGPVIQMYTNVDQDLVGWPAPQGSRSLTPCTCGSSDLYLVTRHADVIPVRRRGDSRGSLLSPRPI"
        "SYLKGSSGGPLLCPAGHAVGLFRAAVCTRGVAKAVDFIPVENLETTMRSPVFTDNSSPPAVPQSFQ"
        "VAHLHAPTGSGKSTKVPAAYAAQGYKVLVLNPSVAATLGFGAYMSKAHGVDPNIRTGVRTITTGSP"
        "ITYSTYGKFLADGGCSGGAYDIIICDECHSTDATSILGIGTVLDQAETAGARLVVLATATPPGSVT"
        "VSHPNIEEVALSTTGEIPFYGKAIPLEVIKGGRHLIFCHSKKKCDELAAKLVALGINAVAYYRGLD"
        "VSVIPTSGDVVVVSTDALMTGFTGDFDSVIDCNTCVTQTVDFSLDPTFTIETTTLPQDAVSRTQRR"
        "GRTGRGKPGIYRFVAPGERPSGMFDSSVLCECYDAGCAWYELTPAETTVRLRAYMNTPGLPVCQDH"
        "LEFWEGVFTGLTHIDAHFLSQTKQSGENFPYLVAYQATVCARAQAPPPSWDQMWKCLIRLKPTLHG"
        "PTPLLYRLGAVQNEVTLTHPITKYIMTCMSADLEVVT"
    )

    def test_ns3_aa_translation(self):
        start, end = GENE_POS['NS3']
        ns3_aa_seq = self.translate_hcv1a(start, end)
        self.assertEqual(
            len(ns3_aa_seq),
            len(self.NS3_AA),
            "Expected reference and translated AA seqs to have equal lengths",
        )
        self.assertEqual(ns3_aa_seq, self.NS3_AA)

    NS5A_AA = (
        "SGSWLRDIWDWICEVLSDFKTWLKAKLMPQLPGIPFVSCQRGYRGVWRGDGIMHTRCHCGAEITGHV"
        "KNGTMRIVGPRTCRNMWSGTFPINAYTTGPCTPLPAPNYKFALWRVSAEEYVEIRRVGDFHYVSGM"
        "TTDNLKCPCQIPSPEFFTELDGVRLHRFAPPCKPLLREEVSFRVGLHEYPVGSQLPCEPEPDVAVL"
        "TSMLTDPSHITAEAAGRRLARGSPPSMASSSASQLSAPSLKATCTANHDSPDAELIEANLLWRQEM"
        "GGNITRVESENKVVILDSFDPLVAEEDEREVSVPAEILRKSRRFARALPVWARPDYNPPLVETWKK"
        "PDYEPPVVHGCPLPPPRSPPVPPPRKKRTVVLTESTLSTALAELATKSFGSSSTSGITGDNTTTSS"
        "EPAPSGCPPDSDVESYSSMPPLEGEPGDPDLSDGSWSTVSSGADTEDVVCC"
    )

    def test_ns5a_aa_translation(self):
        start, end = GENE_POS['NS5A']
        ns5a_aa_seq = self.translate_hcv1a(start, end)
        self.assertEqual(
            len(ns5a_aa_seq),
            len(self.NS5A_AA),
            "Expected reference and translated AA seqs to have equal lengths",
        )


class TestRandomness(unittest.TestCase):

    def test_random_substitution(self):
        for _ in range(1000):
            orig = random.choice(list(NUCLEOTIDES))
            mut = random_substitution(orig)
            self.assertNotEqual(
                orig,
                mut,
                "Expected `random_mutation` to produce a different nucleotide",
            )

    def test_apply_mutation(self):
        seq = functools.reduce(
            lambda seq, _: seq + random_substitution(''),
            range(10000),
            ""
        )
        for _ in range(100):
            gene = random.choice(["NS3", "NS5A", "NS5B"])
            mut = random_mutation(seq, gene=gene)
            applied = apply_mutation(mut, seq)
            self.assertEqual(
                mut.org_nt,
                seq[mut.nt_pos],
                "Expected mutation to contain the original nucleotide",
            )
            self.assertEqual(
                mut.sub_nt,
                applied[mut.nt_pos],
                "Expected new sequence to contain the mutation",
            )

    def test_codon_at(self):
        seq = "abcdef"
        cases = [
            (0, "abc"),
            (1, "abc"),
            (2, "abc"),
            (3, "def"),
            (4, "def"),
            (5, "def"),
        ]
        for idx, expected in cases:
            self.assertEqual(
                codon_at(idx, seq),
                expected,
            )


# ---------------------------------------------------------------------

def make_mut_in_file(infile, outfile='hcv1a_mut.fasta'):
    with open(infile, 'r') as inf:
        inf_header, inf_seq = inf.readlines()
    gene = random.choice(["NS3", "NS5A", "NS5B"])
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

    def check_simple_sub(self):
        # Put the RNG in a known state so we can reproduce test
        # failures
        seed = secrets.token_hex()
        random.seed(seed)

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
            print(seed)
            raise e
        mtn = next(iter(report['mutations']))
        try:
            self.check_vmtn_mtn(vmtn, mtn)
        except AssertionError as e:
            print(e)
            print(mtn)
            print(aln)
            print(seed)
            raise e

    # @unittest.skip("Takes too long")
    def test_simple_subs(self):
        print()
        for i in range(self.iterations):
            if i % 20 == 0 and i > 0:
                print(i)
            self.check_simple_sub()
