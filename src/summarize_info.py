import argparse
import numpy as np
import pysam
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [None]'
    )
    parser.add_argument(
        '-i', '--info', default='AF',
        help='Info tags to query. Tags are separated by commas (e.g. AF,AC,AD). [AF]'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output TSV. ["-"]'
    )
    args = parser.parse_args()
    return args


def summarize_info(fn_vcf, info, fn_out):
    f_vcf = pysam.VariantFile(fn_vcf)
    info = info.split(',')
    if len(info) != len(set(info)):
        raise ValueError('Error: repetitive tags in `info`. Please check.')
        
    info_summ = {}
    info_summ['CHR'] = []
    info_summ['POS'] = []
    info_summ['REF'] = []
    info_summ['ALT'] = []
    for tag in info:
        info_summ[tag] = []
    for var in f_vcf.fetch():
        for i, alt in enumerate(var.alts):
            info_summ['CHR'].append(var.contig)
            info_summ['POS'].append(var.start + 1)
            info_summ['REF'].append(var.alleles[0])
            info_summ['ALT'].append(alt)
            for tag in info:
                info_summ[tag].append(var.info[tag][i])
    df = pd.DataFrame.from_dict(info_summ)
    df.to_csv(fn_out, sep='\t', float_format='%.6f', index=None)
    print(np.histogram(info_summ['AF'], bins=100))


if __name__ == '__main__':
    args = parse_args()
    summarize_info(fn_vcf=args.vcf, info=args.info, fn_out=args.out)