'''
Perform set operations on two VCFs.
E.g. Finding the intersection/union set.
'''

import argparse
import pysam
import sys
from allele_match import fetch_nearby_variants


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [None]'
    )
    parser.add_argument(
        '-p', '--panel',
        help='Path to the reference panel VCF (TBI or CSI indexes are required). [None]'
    )
    parser.add_argument(
        '-r', '--ref',
        help='Path to the reference FASTA (FAI index is required). [None]'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output VCF. Set to "-" to print to stdout. ["-"]'
    )
    parser.add_argument(
        '--happy', action='store_true',
        help='Set for hap.py VCFs. [False]'
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Set to print debug messages. [False]'
    )
    args = parser.parse_args()
    return args



def annotate_vcf(
    fn_vcf: str, fn_panel_vcf: str, fn_fasta: str, fn_out: str, happy_vcf: bool, debug: bool=False
    # af_cutoff: float=0, af_prefix: str=None
    ) -> None:
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        f_vcf.header.add_meta(
            'INFO',
            items=[('ID','AF'), ('Number','A'),
                   ('Type','Float'), ('Description','Population allele frequency')])
    except:
        raise ValueError(f'Error: Cannot open "{fn_vcf}"')
    try:
        f_panel = pysam.VariantFile(fn_panel_vcf)
    except:
        raise ValueError(f'Error: Cannot open "{fn_panel_vcf}"')
    try:
        f_fasta = pysam.FastaFile(fn_fasta)
    except:
        raise ValueError(f'Error: Cannot open "{fn_fasta}"')
    try:
        f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
    except:
        raise ValueError(f'Error: Cannot create "{fn_out}"')

    # if af_cutoff > 0 and af_prefix == None:
    #     raise ValueError(f'Error: `allele-frequency-prefix` needs to be set when `allele-frequency-cutoff` > 0')
    # elif af_cutoff > 0:
    #     f_out_high = pysam.VariantFile(af_prefix+f'-af_gt_{af_cutoff}.vcf', 'w', header=f_vcf.header)
    #     f_out_low = pysam.VariantFile(af_prefix+f'-af_leq_{af_cutoff}.vcf', 'w', header=f_vcf.header)
    
    for var in f_vcf.fetch():
        if happy_vcf:
            # Only check variants in confident regions (hap.py specific)
            if var.info.get('Regions'):
                fetch_nearby_cohort(var, f_panel, f_fasta, f_out, debug)
        else:
            if var.filter.get('PASS'):
                # Only take 'PASS' variants
                fetch_nearby_cohort(var, f_panel, f_fasta, f_out, debug)


if __name__ == '__main__':
    # We use python3.6 because we require dict is ordered.
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()

    annotate_vcf(
        fn_vcf=args.vcf,
        fn_panel_vcf=args.panel,
        fn_fasta=args.ref,
        fn_out=args.out,
        happy_vcf=args.happy,
        debug=args.debug
        # af_cutoff=args.allele_frequency_cutoff,
        # af_prefix=args.allele_frequency_prefix
    )
