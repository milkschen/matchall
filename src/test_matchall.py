import matchall
import filter
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


class TestMatchAll(unittest.TestCase):
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
        var = matchall.match_allele(
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
    def test_match_matchall(
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
        var = matchall.match_allele(
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
            # print(v.samples.keys())
            # print(v.samples[v.samples.keys()[0]])
            v_result = matchall.fetch_nearby_cohort(
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
    
    def test_vmerge_variant_pair_no_gt(self):
        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add('chr1')
        # info = {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Population allele frequency'}
        # vcf_header.add_meta('INFO', items=info.items())
        # Make the variant record under test
        v1 = vcf_header.new_record(
            contig='chr1',
            start=10814,
            stop=10815,
            alleles=('T', 'TC'))
        v2 = vcf_header.new_record(
            contig='chr1',
            start=10815,
            stop=10816,
            alleles=('C', 'CCA'))
        v = matchall.vmerge_variant_pair(v1=v1, v2=v2)
        self.assertEqual(v.start, 10814)
        self.assertEqual(v.stop, 10816)
        self.assertEqual(v.alleles, ('T', 'TCCA'))

    def test_vmerge_variant_pair_diploid(self):
        # TODO
        f_vcf = pysam.VariantFile(os.path.join('test_data/tmp', 'HG003.octopus.overlapped.chr1.vcf'))
        for i, v in enumerate(f_vcf):
            if i == 0:
                v1 = v
            elif i == 1:
                v2 = v
            else:
                break
        v = matchall.vmerge_variant_pair(v1=v1, v2=v2)
        self.assertEqual(v.start, 10814)
        self.assertEqual(v.stop, 10816)
        self.assertEqual(v.alleles, ('T', 'TC', 'TCCA'))


class TestFilter(unittest.TestCase):
    def setUp(self):
        pass

    @parameterized.expand([
        [
            'AF > 0.5',
            ('>', ['AF', 0.5], None)
        ],
        [
            'GT[0] = 0|1',
            ('=', ['GT', "0|1"], 0)
        ],
        [
            'BD!=FP',
            ('!=', ['BD', 'FP'], None)
        ]
    ])
    def test_filter_parse_expression(self, exp, gold):
        result = filter.parse_expression(exp)
        opt = result[0]
        comps = result[1]
        index = result[2]
        self.assertEqual(opt, gold[0])
        self.assertListEqual(comps, gold[1])
        self.assertEqual(index, gold[2])

    @parameterized.expand([
        ['>', 1, 2, False],
        ['!=', 5, 4, True],
        ['==', 3.1, 3.1, True],
        ['=', 'FP', 'FP', True]
    ])
    def test_eval_operator(self, operator, left, right, gold):
        result = filter.eval_operator(operator, left, right)
        self.assertEqual(result, gold)

    @parameterized.expand([
        [
            PysamVariant(
                start=100, stop=101, alleles=('C', 'A')
            ),
            ('AF', (0.09824279695749283,)),
            'AF<0.1',
            True
        ],
        [
            PysamVariant(
                start=100, stop=101, alleles=('C', 'A')
            ),
            ('AF', (0.09824279695749283,)),
            'AF>=0.1',
            False
        ]
    ])
    def test_vcf_filter_core_info(
        self, var, info_pattern, info_exps, gold):
        vcf_header = pysam.VariantHeader()
        vcf_header.contigs.add('chr1')
        info = {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Population allele frequency'}
        # vcf_header.formats.add('BD', '1', 'String', 'BD')
        vcf_header.add_meta('INFO', items=info.items())
        record = vcf_header.new_record(
            contig=var.contig,
            start=var.start,
            stop=var.stop,
            alleles=var.alleles)
        record.info.__setitem__(info_pattern[0], info_pattern[1])
        ie = None
        if info_exps:
            ie = [filter.parse_expression(info_exps)]
        result = filter.vcf_filter_core(
            var=record, info_exps=ie, format_exps=None
        )
    
    @parameterized.expand([
        ["DP>35", [True, True, True, False, True, True, True, True, False, False]],
        ["GQ>20", [True, True, True, True, True, False, True, True, True, False]]
    ])
    def test_vcf_filter_core_format(self, format_exp, gold):
        f_test = pysam.VariantFile(
            os.path.join('test_data', 'HG00733-hifi_deepvariant-chr20_568936_571052.vcf.gz'))
        results = []
        fe = [filter.parse_expression(format_exp)]
        for var in f_test.fetch():
            results.append(filter.vcf_filter_core(var=var, info_exps=None, format_exps=fe))
        self.assertListEqual(results, gold)


if __name__ == '__main__':
    unittest.main()
