# import matchall
import os
import subprocess
import unittest

from parameterized import parameterized
from test_matchall import PysamVariant



class TestAnnotate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fn_vcf = os.path.join('test_data', 'HG00733-hifi_deepvariant-chr20_568936_571052.vcf.gz')
        cls.fn_query_vcf = os.path.join('test_data', 'chr20_560000_580000.cohort.vcf.gz')
        cls.fn_fasta = os.path.join('test_data', 'chr20_1_580000.fa')
        cls.fn_gold = {
            'AF': os.path.join('test_data', 'gold-HG00733-hifi_deepvariant-chr20_568936_571052.af.vcf'),
            'HWE': os.path.join('test_data', 'gold-HG00733-hifi_deepvariant-chr20_568936_571052.hwe.vcf'),
            'NS': os.path.join('test_data', 'gold-HG00733-hifi_deepvariant-chr20_568936_571052.ns.vcf'),
        }

    @parameterized.expand([
        ['AF'],
        ['HWE'],
        ['NS']
    ])
    def test_annotate(self, info_id):
        # Using grep to exclude VCF headers
        cmd = (f'python src/annotate.py -r {self.fn_fasta} -v {self.fn_vcf} '
               f'-q {self.fn_query_vcf} --info {info_id} | grep "^[^#]"')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        gold = b''
        with open(self.fn_gold[info_id], 'rb') as f_gold:
            for line in f_gold:
                gold += line
        self.assertEqual(stdout, gold)


class TestMatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fn_vcf = os.path.join('test_data', 'HG002-chr1_870000_875000.vcf.gz')
        cls.fn_query_vcf = os.path.join('test_data', 'HG002-major-chr1_870000_875000.vcf.gz')
        cls.fn_fasta = os.path.join('test_data', 'chr1_1_880000.fa')
        cls.fn_gold = os.path.join('test_data', 'gold-HG002-chr1_870000_875000.vcf')

    def test_match(self):
        # Using grep to exclude VCF headers
        cmd = (f'python src/compare.py -r {self.fn_fasta} -v {self.fn_vcf} '
               f'-q {self.fn_query_vcf} | grep "^[^#]"')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        gold = b''
        with open(self.fn_gold, 'rb') as f_gold:
            for line in f_gold:
                gold += line
        self.assertEqual(stdout, gold)


class TestIsecPrivate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fn_vcf = os.path.join('test_data', 'chr20_156877-A.vcf.gz')
        cls.fn_query_vcf = os.path.join('test_data', 'chr20_156877-B.vcf.gz')
        cls.fn_fasta = os.path.join('test_data', 'chr20_1_580000.fa')
        cls.fn_gold_isec0 = os.path.join('test_data', 'chr20_156877-isecA.vcf.gz')
        cls.fn_gold_isec1 = os.path.join('test_data', 'chr20_156877-isecB.vcf.gz')
        cls.fn_gold_private0 = os.path.join('test_data', 'chr20_156877-privateA.vcf.gz')
        cls.fn_gold_private1 = os.path.join('test_data', 'chr20_156877-privateB.vcf.gz')
        cls.output_label = 'test-compare-chr20_156877'
        cls.fn_results_isec0 = os.path.join('test_data', f'{cls.output_label}.isec0.vcf.gz')
        cls.fn_results_isec1 = os.path.join('test_data', f'{cls.output_label}.isec1.vcf.gz')
        cls.fn_results_private0 = os.path.join('test_data', f'{cls.output_label}.private0.vcf.gz')
        cls.fn_results_private1 = os.path.join('test_data', f'{cls.output_label}.private1.vcf.gz')
        cmd = (f'python src/compare.py -r {cls.fn_fasta} -v {cls.fn_vcf} '
               f'-q {cls.fn_query_vcf} -m isec,private -op test_data/{cls.output_label}')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
 
    @classmethod
    def tearDownClass(cls):
        cmd = (f'rm {cls.fn_results_isec0};'
               f'rm {cls.fn_results_isec1};'
               f'rm {cls.fn_results_private0};'
               f'rm {cls.fn_results_private1};')
        process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process_rm.communicate()

    def test_isec0(self):
        gold = b''
        with open(self.fn_gold_isec0, 'rb') as f_gold_isec0:
            for line in f_gold_isec0:
                gold += line
        result = b''
        with open(self.fn_results_isec0, 'rb') as f_result_isec0:
            for line in f_result_isec0:
                result += line
        self.assertEqual(result, gold)
    
    def test_isec1(self):
        gold = b''
        with open(self.fn_gold_isec1, 'rb') as f_gold_isec1:
            for line in f_gold_isec1:
                gold += line
        result = b''
        with open(self.fn_results_isec1, 'rb') as f_result_isec1:
            for line in f_result_isec1:
                result += line
        self.assertEqual(result, gold)
    
    def test_private0(self):
        gold = b''
        with open(self.fn_gold_private0, 'rb') as f_gold_private0:
            for line in f_gold_private0:
                gold += line
        result = b''
        with open(self.fn_results_private0, 'rb') as f_result_private0:
            for line in f_result_private0:
                result += line
        self.assertEqual(result, gold)

    def test_private1(self):
        gold = b''
        with open(self.fn_gold_private1, 'rb') as f_gold_private1:
            for line in f_gold_private1:
                gold += line
        result = b''
        with open(self.fn_results_private1, 'rb') as f_result_private1:
            for line in f_result_private1:
                result += line
        self.assertEqual(result, gold)


if __name__ == '__main__':
    unittest.main()