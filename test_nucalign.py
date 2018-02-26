import collections
import functools
import json
import subprocess
import random
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

    def __init__(self, infilename):
        self.nuc_result = self.nucalign(infilename)

    def mutations(self, gene='NS3'):
        return [result['Report']['Mutations'] for result in self.nuc_result[gene]]


class TestReferenceSequence(unittest.TestCase):

    reference_filename = 'hcv1a.fasta'

    def test_no_mutations_in_reference(self):
        alignment = NucAlignment(self.reference_filename)
        muts = alignment.mutations()
        for mut in muts:
            self.assertEqual(len(mut), 0)


# ---------------------------------------------------------------------

GENE_POS = {
    'NS3': (3419, 5312),
    'NS5A': (6257, 7601),
    'NS5B': (7602, 9374),
}

NUCLEOTIDES = {'G', 'C', 'T', 'A'}


Mutation = collections.namedtuple(
    "Mutation",
    ["pos", "sub_nt", "org_nt", "sub_cod", "org_cod"],
)


def random_gene_pos(gene='NS3'):
    return random.choice(range(*GENE_POS[gene]))


def random_substitution(orig):
    return random.choice(list(NUCLEOTIDES - {orig}))


def codon_idx(idx):
    return idx - (idx % 3)


def codon_at(idx, seq):
    cod_idx = codon_idx(idx)
    return seq[cod_idx:cod_idx + 3]


def random_mutation(seq, gene='NS3'):
    idx = random_gene_pos(gene)
    orig_nt = seq[idx]
    orig_codon = codon_at(idx, seq)
    sub = random_substitution(orig_nt)
    sub_codon = orig_codon[:idx%3] + sub + orig_codon[(idx%3)+1:]
    return Mutation(
        pos=idx,
        org_nt=orig_nt,
        org_cod=orig_codon,
        sub_nt=sub,
        sub_cod=sub_codon,
    )


def apply_mutation(mut, sequence):
    idx = mut.pos
    return (
        sequence[:idx]
        + mut.sub_nt
        + sequence[idx+1:]
    )


class TestTranslation(unittest.TestCase):

    @staticmethod
    def translate_hcv1a(start, end):
        _, hcv1a_nt_seq = open("hcv1a.fasta").readlines()
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

    def test_gene_pos(self):
        for _ in range(1000):
            ridx = random_gene_pos('NS3')
            self.assertIn(ridx, range(3420, 5474))

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
            mut = random_mutation(seq)
            applied = apply_mutation(mut, seq)
            self.assertEqual(
                mut.org_nt,
                seq[mut.pos],
                "Expected mutation to contain the original nucleotide",
            )
            self.assertEqual(
                mut.sub_nt,
                applied[mut.pos],
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
    mutation = random_mutation(inf_seq)
    idx = mutation.pos
    outf_header = inf_header.strip() + " with mutation: {}\n".format(mutation)
    outf_seq = apply_mutation(mutation, inf_seq)
    with open(outfile, 'w') as outf:
        outf.writelines([outf_header, outf_seq])
    return mutation


class TestSimpleSubstitutions(unittest.TestCase):

    reference_file = "hcv1a.fasta"
    tmp_filename = "hcv1a_mut.fasta"
    iterations = 10

    def check_simple_sub(self):
        mutation = make_mut_in_file(
            self.reference_file,
            outfile=self.tmp_filename,
        )
        alignment = NucAlignment(self.tmp_filename)
        return mutation, alignment

    def test_simple_sub(self):
        print()
        hdr = "VirtMut         FoundMut            Diff"
        print(hdr)

        def format_aln_mutation(m):
            if m is None:
                return ""
            return "{} {} {} (nt:{})".format(
                m["ReferenceText"],
                m["Position"],
                m["AminoAcidText"],
                m["NAPosition"],
            )

        def format_virt_mutation(mtn):
            return "{} {} {}".format(
                mtn.org_nt.lower(),
                mtn.pos,
                mtn.sub_nt.lower(),
            )

        def pos_diff(virt, aln):
            if virt is None or aln is None:
                return "-"
            return virt.pos - aln["NAPosition"]

        def format_codons(org, sub):
            return "{}({}) -> {}({})".format(
                org,
                Bio.Seq.translate(org),
                sub,
                Bio.Seq.translate(sub),
            )

        tmpl = "{virtmut: <16}{foundmut: <20}{diff: <16}"
        for i in range(1):
            mtn, aln = self.check_simple_sub()
            mutations = aln.mutations()
            assert len(mutations) == 1
            mutations = mutations[0]  # 😢
            assert len(mutations) < 2
            mutation = next(iter(mutations), None)
            data = {
                "virtmut": format_virt_mutation(mtn),
                "foundmut": format_aln_mutation(mutation),
                "diff": pos_diff(mtn, mutation),
            }
            outp = tmpl.format(**data)
            print(outp)
