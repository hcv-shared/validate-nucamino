import collections
import functools
import json
import os
import random
import secrets
import subprocess
import unittest

import Bio.Seq


def print_seed_on_assertionerror(f):

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        seed = secrets.token_hex()
        random.seed(seed)
        try:
            f(*args, **kwargs)
        except AssertionError as e:
            print("seed = {}".format(seed))
            raise e

    return wrapped


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
# Mutations (SNP)

GENES = ["NS3", "NS5A", "NS5B"]
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


def random_mutation(seq, gene=None):
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
            gene = random.choice(GENES)
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

    @unittest.skip("Takes too long")
    def test_simple_subs(self):
        print()
        for i in range(self.iterations):
            if i % 20 == 0 and i > 0:
                print(i)
            self.check_simple_sub()


# ---------------------------------------------------------------------
# Generate and Work with Indels

def random_nt():
    return random.choice("GATC")


def apply_insertion(seq, pos, ins):
    return seq[:pos] + ins + seq[pos:]


def apply_deletion(seq, pos, count):
    return seq[:pos] + seq[pos+count:]


def random_insertion(seq, max_length=1, gene=None):
    start, end = GENE_POS[gene]
    ins_length = random.randint(1, max_length)
    pos = start + random.randint(1, (end - start - ins_length))
    ins = "".join(random_nt() for _ in range(ins_length))
    return (pos, ins)


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
        for pos, ins, expected in cases:
            inserted = apply_insertion(base_seq, pos, ins)
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
        start, end = GENE_POS[gene]
        rand_pos, rand_ins = random_insertion(
            "",
            max_length=max_len,
            gene=gene,
        )
        self.assertLessEqual(start, rand_pos)
        self.assertLessEqual(rand_pos, end)
        self.assertLessEqual(len(rand_ins), max_len)
        self.assertTrue(set(rand_ins).issubset(NUCLEOTIDES))

    @print_seed_on_assertionerror
    def test_insertion_generation(self):
        for i in range(100000):
            self.check_insertion_generation()

    # TODO(nknight): test deltion generation

# ---------------------------------------------------------------------
# Test Alinger on Indels

class TestInsertions(unittest.TestCase):


    def setUp(self):
        with open("hcv1a.fasta") as inf:
            _, self.hcv1a_seq = inf.readlines()

    def create_insertion_file(self):
        gene = 'NS3'
        ins = random_insertion(
            self.hcv1a_seq,
            max_length=1,
            gene=gene,
        )
        applied = apply_insertion(self.hcv1a_seq, *ins)
        with open("hcv1a_mut.fasta", "w") as outf:
            outf.writelines([
                "> HCV1A with insertion {} {}\n".format(*ins),
                applied,
            ])

    @print_seed_on_assertionerror
    def test_make_insertion(self):
        self.create_insertion_file()
