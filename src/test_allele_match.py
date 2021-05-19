import allele_match
import os
import pysam
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
            'bi-allelic var, one cohort',
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
        [ 
            'bi-allelic var, two cohorts',
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
        [ 
            'tri-allelic var, multiple cohorts',
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
        [ 
            'bi-allelic var, a complex cohort',
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
        ],
        [ 
            'a special case I only see from the octopus variant caller',
            PysamVariant(
                start=143327539, stop=143327543, alleles=('CACA', '*', '*')
            ),
            [
                PysamVariant(
                    start=143327528, stop=143327551, alleles=('GACAGCGGCGGCACAGCGGCGGC', 'GACAGCGGCGGC', 'GACAGCGGCGGCGCGGCGGC', 'GATAGCGGCGGCACAGCGGCGGC', 'GACA', 'GCCAGCGGCGGCACAGCGGCGGC'),
                    af=(0.0692891,0.000998403,0.00519169,0.000199681,0.000199681)
                )
            ],
            'GACAGCGGCGGCACAGCGGCGGC',
            (0,0)
        ]
    ])
    def test_match_allele_af(
        self, description, var, cohorts, ref, gold):
        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add('chr1')
        info = {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Population allele frequency'}
        vcf_header.add_meta('INFO', items=info.items())
        # Make the variant record under test
        record = vcf_header.new_record(
            contig=var.contig,
            start=var.start,
            stop=var.stop,
            alleles=var.alleles)
        # Make cohort variants
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
        ref = ref
        var = allele_match.match_allele(
            var=record, cohort_vars=c_records,
            ref=ref, update_info=info, query_info=info)
        
        for i, a in enumerate(var.info['AF']):
            self.assertAlmostEqual(a, gold[i])

    @parameterized.expand([
        [ 
            'bi-allelic var, matched',
            PysamVariant(
                start=19003, stop=19005, alleles=('AA', 'AG')
            ),
            [
                PysamVariant(
                    start=19004, stop=19005, alleles=('A', 'G')
                )
            ],
            'AA',
            (1,)
        ],
        [   
            'bi-allelic var, unmatched',
            PysamVariant(
                start=19003, stop=19005, alleles=('AA', 'AG')
            ),
            [
                PysamVariant(
                    start=19003, stop=19004, alleles=('A', 'G')
                )
            ],
            'AA',
            (0,)
        ],
        [ 
            'one allele is matched in a multi-allelic site',
            PysamVariant(
                start=6540123, stop=6540132, alleles=('C', 'CTTTTT', 'CTTTTTTTT')
            ),
            [
                PysamVariant(
                    start=6540123, stop=6540132, alleles=('C', 'CTTTTTTTT')
                )
            ],
            'CTTTTTTTT',
            (0, 1)
        ]
    ])
    def test_match_allele_match(
        self, description, var, cohorts, ref, gold):
        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add('chr1')
        info = {'ID': 'MATCH', 'Number': 'A', 'Type': 'Integer', 'Description': 'If genotype is matched with a query'}
        # info = {'ID': 'MATCH', 'Number': '1', 'Type': 'Integer', 'Description': 'If genotype is matched with a query'}
        vcf_header.add_meta('INFO', items=info.items())
        # Make the variant record under test
        record = vcf_header.new_record(
            contig=var.contig,
            start=var.start,
            stop=var.stop,
            alleles=var.alleles)
        # Make cohort variants
        c_records = []
        for cohort in cohorts:
            v = vcf_header.new_record(
                contig = cohort.contig,
                start = cohort.start,
                stop = cohort.stop,
                alleles = cohort.alleles
            )
            # v.info.__setitem__('AF', cohort.af)
            c_records.append(v)
        ref = ref
        var = allele_match.match_allele(
            var=record, cohort_vars=c_records,
            ref=ref, update_info=info, query_info=None)
        
        for i, a in enumerate(var.info['MATCH']):
            self.assertAlmostEqual(a, gold[i])
        # self.assertAlmostEqual(var.info['MATCH'], gold)

    @parameterized.expand([
        [ # AF, Number='A'
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'},
            [0.009784, 0.539137, 0.170927, 0.231030, 0.231030, 
             0.184105, 0.224641, 0.485423, 0.539337, 0.243211]
        ],
        [ # NS, Number='1'
            {'ID': 'NS', 'Number': '1', 'Type': 'Integer', 'Description': 'Number of samples with data'},
            [2504 for _ in range(10)]
        ]
    ])
    def test_fetch_nearby_cohort(self, info, gold):
        f_vcf = pysam.VariantFile(os.path.join('test_data', 'HG00733-hifi_deepvariant-chr20_568936_571052.vcf.gz'))
        f_query_vcf = pysam.VariantFile(os.path.join('test_data', 'chr20_560000_580000.cohort.vcf.gz'))
        f_fasta = pysam.FastaFile(os.path.join('test_data', 'chr20_1_580000.fa'))
        f_vcf.header.add_meta('INFO', items=info.items())

        for i, v in enumerate(f_vcf.fetch()):
            v_result = allele_match.fetch_nearby_cohort(
                var=v, f_query_vcf=f_query_vcf,
                f_fasta=f_fasta, update_info=info, query_info=info)
            
            # Number of values = number of alt alleles
            # Here all tuples are unit-length
            if info['Number'] == 'A':
                result = v_result.info[info['ID']][0]
            # Number of values = 1
            elif info['Number'] == '1':
                result = v_result.info[info['ID']]

            self.assertAlmostEqual(result, gold[i], places=6)


if __name__ == '__main__':
    unittest.main()