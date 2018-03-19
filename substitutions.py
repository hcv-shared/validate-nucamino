import collections
import functools
import random
import unittest

import Bio.Seq

import common

_Mutation = collections.namedtuple(
    "_Mutation",
    [
        "nt_pos",
        "aa_pos",
        "sub_nt",
        "org_nt",
        "sub_cod",
        "org_cod",
        "gene",
        "genotype",
    ],
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


def gene_seq(seq, genotype, gene):
    start, end = common.GENE_POS[genotype][gene]
    return seq[start:end]


def random_substitution(orig):
    return random.choice(list(common.NUCLEOTIDES - {orig}))


def codon_idx(idx):
    return idx - (idx % 3)


def codon_at(idx, seq):
    cod_idx = codon_idx(idx)
    return seq[cod_idx:cod_idx + 3]


def random_mutation(seq, gene, genotype):
    start, end = common.GENE_POS[genotype][gene]
    gseq = gene_seq(seq, gene=gene, genotype=genotype)
    idx = random.choice(range(len(gseq)))
    orig_nt = gseq[idx]
    orig_codon = codon_at(idx, gseq)
    sub = random_substitution(orig_nt)
    sub_codon = orig_codon[:idx % 3] + sub + orig_codon[(idx % 3)+1:]
    return Mutation(
        nt_pos=idx+start,
        aa_pos=idx//3 + 1,
        org_nt=orig_nt,
        org_cod=orig_codon,
        sub_nt=sub,
        sub_cod=sub_codon,
        gene=gene,
        genotype=genotype,
    )


def apply_mutation(mut, sequence):
    idx = mut.nt_pos
    return (
        sequence[:idx]
        + mut.sub_nt
        + sequence[idx+1:]
    )


# ---------------------------------------------------------------------
# Tests

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
        start, end = common.GENE_POS['1a']['NS3']
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
        start, end = common.GENE_POS['1a']['NS5A']
        ns5a_aa_seq = self.translate_hcv1a(start, end)
        self.assertEqual(
            len(ns5a_aa_seq),
            len(self.NS5A_AA),
            "Expected reference and translated AA seqs to have equal lengths",
        )


class TestRandomness(unittest.TestCase):

    def test_random_substitution(self):
        for _ in range(1000):
            orig = random.choice(list(common.NUCLEOTIDES))
            mut = random_substitution(orig)
            self.assertNotEqual(
                orig,
                mut,
                "Expected random_substitution to yield a different nucleotide",
            )

    def test_apply_mutation(self):
        seq = functools.reduce(
            lambda seq, _: seq + random_substitution(''),
            range(10000),
            ""
        )
        for _ in range(100):
            gene = random.choice(common.GENES)
            genotype = random.choice(common.GENOTYPES)
            mut = random_mutation(seq, gene=gene, genotype=genotype)
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
