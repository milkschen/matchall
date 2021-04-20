import allele_match
import pysam
# import pytest
import unittest
from parameterized import parameterized

# @pytest.fixture
# def vcf_header():
#     vcf_header = pysam.VariantHeader()
#     vcf_header.add_sample("sample1")
#     # vcf_header.add_sample("sample2")
#     vcf_header.contigs.add("chr1")
#     return vcf_header

class PysamVariant():
    def __init__(self, start, stop, alleles, contig='chr1', af=None):
        self.contig = contig
        self.start = start
        self.stop = stop
        self.alleles = alleles
        self.af = af


class TestAlleleMatch(unittest.TestCase):
    def setUp(self):
        pass
    
    @parameterized.expand([
        [ # bi-allelic var, one cohort
            PysamVariant(
                start=19003, stop=19004, alleles=('A', 'G')
            ),
            [
                PysamVariant(
                    start=19003, stop=19004, alleles=('A', 'G'),
                    af=0.401558
                )
            ],
            'A',
            (0.401558,)
        ],
        [ # bi-allelic var, two cohorts
            PysamVariant(
                start=20093, stop=20096, alleles=('TAA', 'T')
            ),
            [
                PysamVariant(
                    start=20093, stop=20098, alleles=('TAAAC', 'T'),
                    af=0.0107827
                ),
                PysamVariant(
                    start=20093, stop=20096, alleles=('TAA', 'T'),
                    af=0.0455272
                )
            ],
            'TAC',
            (0.0455272,)
        ],
        [ # tri-allelic var, multiple cohorts
            PysamVariant(
                start=181582, stop=181587, alleles=('CGGGG', 'C', 'CG')
            ),
            [
                PysamVariant(
                    start=181582, stop=181589, alleles=('CGGGGGG', 'CGG', 'CGGG', 'CGGGGG', 'CGGGG', 'CG', 'CGGGGGGG', 'CGGGGGGGG', 'C', 'CGGGGGGGGG', 'CGGGGGC', 'CGGGGAGG'),
                    af=(0.333866,0.113818,0.0810703,0.0513179,0.0171725,0.013778,0.00259585,0.00259585,0.00119808,0.000199681,0.000199681)
                ),
                PysamVariant(
                    start=181582, stop=181591, alleles=('CGGGGGGGG', 'C'),
                    af=0.00539137
                )
            ],
            'CGGGGGGGGG',
            (0.333866, 0.113818)
        ],
        [ # bi-allelic var, a complex cohort
            PysamVariant(
                start=893788, stop=893825, alleles=('AAAAAAAAAAAAAATATATATATATATATATATATAT', 'A'),
            ),
            [
                PysamVariant(
                    start=893786, stop=893827, alleles=('AAAAAAAAAAAAAAAATATATATATATATATATATATATAT', 'AAAAT', 'AATAT', 'AAAAA', 'AAAAAAAAAAAAAATATATATATATATATATATATATATAT', 'AAAAAAAAAAAAAATATATAT', 'AAAAAAAAAAAAAAAAAATATATATATATATATATATATAT', 'AAAAAAAAAAAAAAAATATATATATATATATAT', 'AAAAAAAAAAAATATATATAT', 'AAAAAAAAAAAAAAAATATATAT', 'AAAAAAAAAAAAAAAATATAT', 'AAAAAAAAAAAAAAATATATATATATATATATATATATATATAT', 'AAAAAAAAAAAAAAAATATATATAT', 'AAAAAAAAAAAAAAAAATATAT', 'AAAAAAAAAAAAAATATATATAT', 'AAAAAAAAAAATATATATATATATATATATATAT', 'AAAAAAAAAAAAAAATATATATATAT', 'AAAAAAAAAAAAATATATATATATAT'),
                    af=(0.734225,0.0365415,0.0361422,0.000998403,0.00219649,0.000998403,0.00119808,0.000998403,0.000798722,0.000599042,0.000399361,0.000199681,0.000199681,0.000199681,0.000199681,0.000199681,0.000199681)
                )
            ],
            'AAAAAAAAAAAAAAAATATATATATATATATATATATATATATA',
            (0.734225,)
        ]
    ])
    def test_compare_haplotype(self, var, cohorts, ref, gold):
        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add('chr1')
        vcf_header.add_meta('INFO', items=[('ID','AF'), ('Number','A'), ('Type','Float'), ('Description','Population allele frequency')])
        record = vcf_header.new_record(
            contig=var.contig,
            start=var.start,
            stop=var.stop,
            alleles=var.alleles)
        c_records = []
        for cohort in cohorts:
            v = vcf_header.new_record(
                contig = cohort.contig,
                start = cohort.start,
                stop = cohort.stop,
                alleles = cohort.alleles
            )
            v.info.__setitem__('AF', cohort.af)
            c_records.append(v)
        ref=ref
        var = allele_match.compare_haplotypes(record, c_records, ref)
        for i, a in enumerate(var.info['AF']):
            self.assertAlmostEqual(a, gold[i])
        # self.assertAlmostEqual(var.info['AF'][0], gold)

if __name__ == '__main__':
    unittest.main()
