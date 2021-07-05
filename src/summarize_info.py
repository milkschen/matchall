''' 
Summarize an info field in a VCF as a TSV. 

Columns:
CHR POS REF ALT TYPE INFO1 INFO2 ... FORMAT1 FORMAT2 ...
'''

import argparse
# import numpy as np
import pandas as pd
import pysam


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [None]'
    )
    parser.add_argument(
        '-i', '--info', default='AF',
        help='Info tags to retrieve. Tags are separated by commas (e.g. AF,AC,AD). [AF]'
    )
    parser.add_argument(
        '-f', '--format', default=None,
        help='Format tags to retrieve. Tags are separated by commas (e.g. BD,GT). [None]'
    )
    parser.add_argument(
        '--format-idx', default=1, type=int,
        help='Index of samples to retrieve format tags. [1]'
    )
    # parser.add_argument(
    #     '-t', '--type', default='all',
    #     help='Variant type. Options: all, snp, indel [all]'
    # )
    parser.add_argument(
        '--happy', action='store_true',
        help='Set for hap.py VCFs. Will only consider variants in confident regions. [False]'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output TSV. ["-"]'
    )
    args = parser.parse_args()
    return args


def summarize_info(
    fn_vcf: str,
    fn_out: str, 
    info: str,
    vcf_format: str,
    format_idx: int,
    happy_vcf: bool
    ) -> None:
    f_vcf = pysam.VariantFile(fn_vcf)
    if info:
        info = info.split(',')
    if vcf_format:
        vcf_format = vcf_format.split(',')
    if len(info) != len(set(info)):
        raise ValueError('Error: repetitive tags in `info`. Please check.')
        
    info_summ = {}
    info_summ['CHR'] = []
    info_summ['POS'] = []
    info_summ['REF'] = []
    info_summ['ALT'] = []
    info_summ['TYPE'] = []
    if info:
        for tag in info:
            info_summ[tag] = []
    if vcf_format:
        for tag in vcf_format:
            info_summ[tag] = []

    for var in f_vcf.fetch():
        if happy_vcf and not var.info.get('Regions'):
            pass
        for i, alt in enumerate(var.alts):
            # Skip when see an empty ALT
            if alt == '*':
                continue
            info_summ['CHR'].append(var.contig)
            info_summ['POS'].append(var.start + 1)
            info_summ['REF'].append(var.alleles[0])
            info_summ['ALT'].append(alt)
            if len(var.alleles[0]) == 1 and len(alt) == 1:
                info_summ['TYPE'].append('SNP')
            else:
                info_summ['TYPE'].append('INDEL')
            if info:
                for tag in info:
                    info_summ[tag].append(var.info[tag][i])
            if vcf_format:
                for tag in vcf_format:
                    value = var.samples[format_idx][tag]
                    if type(value) is tuple:
                        info_summ[tag].append(value[i])
                    else:
                        info_summ[tag].append(value)

    df = pd.DataFrame.from_dict(info_summ)
    df.to_csv(fn_out, sep='\t', float_format='%.6f', index=None)
    # print(np.histogram(info_summ['AF'], bins=100))


if __name__ == '__main__':
    args = parse_args()
    summarize_info(
        fn_vcf=args.vcf, fn_out=args.out, 
        info=args.info, vcf_format=args.format, 
        format_idx=args.format_idx, happy_vcf=args.happy
        )
