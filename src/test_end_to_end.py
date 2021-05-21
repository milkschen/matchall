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
        cmd = (f'python src/set_operation.py -r {self.fn_fasta} -v {self.fn_vcf} '
               f'-q {self.fn_query_vcf} | grep "^[^#]"')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        gold = b''
        with open(self.fn_gold, 'rb') as f_gold:
            for line in f_gold:
                gold += line
        self.assertEqual(stdout, gold)


if __name__ == '__main__':
    unittest.main()