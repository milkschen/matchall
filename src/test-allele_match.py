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
        [
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
            'TAC'
        ]
    ])
    def test_compare_haplotype(self, var, cohorts, ref):
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
        # record = vcf_header.new_record(
        #     contig='chr1',
        #     start=20093,
        #     stop=20096,
        #     alleles=('TAA', 'T'))
        # cohorts = [
        #     vcf_header.new_record(
        #         # contig='chr1',
        #         start=20093,
        #         stop=20098,
        #         alleles=('TAAAC', 'T')),
        #     vcf_header.new_record(
        #         # contig='chr1',
        #         start=20093,
        #         stop=20096,
        #         alleles=('TAA', 'T'))
        # ]
        # cohorts[0].info.__setitem__('AF', 0.0107827)
        # cohorts[1].info.__setitem__('AF', 0.0455272)
        # ref='TAC'
        var = allele_match.compare_haplotypes(record, c_records, ref)
        # print(var.info['AF'][0])
        self.assertAlmostEqual(var.info['AF'][0], 0.0455272)

if __name__ == '__main__':
    unittest.main()
