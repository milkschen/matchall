'''
Perform set operations on two VCFs.
E.g. Finding the intersection/union set.
'''
import allele_match
import argparse
import pysam
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf', required=True,
        help='Path to target VCF. [required]'
    )
    parser.add_argument(
        '-q', '--query-vcf', required=True,
        help='Path to query VCF (TBI or CSI indexes are required). [required]'
    )
    parser.add_argument(
        '-r', '--ref', required=True,
        help='Path to the reference FASTA (FAI index is required). [required]'
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


def match_vcf(
    fn_vcf: str, fn_query_vcf: str, fn_fasta: str,
    fn_out: str, happy_vcf: bool, debug: bool=False
    ) -> None:
    update_info = {'ID': 'MATCH', 'Number': '1', 'Type': 'Integer', 'Description': 'If genotype is matched with a query'}
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        f_vcf.header.add_meta(
            'INFO',
            items=update_info.items())
    except:
        raise ValueError(f'Error: Cannot open "{fn_vcf}"')

    try:
        f_query_vcf = pysam.VariantFile(fn_query_vcf)
    except:
        raise ValueError(f'Error: Cannot open "{fn_query_vcf}"')

    try:
        f_fasta = pysam.FastaFile(fn_fasta)
    except:
        raise ValueError(f'Error: Cannot open "{fn_fasta}"')
        
    try:
        f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
    except:
        raise ValueError(f'Error: Cannot create "{fn_out}"')

    def select_variant(var):
        if happy_vcf and var.info.get('Regions'):
            # Check variants in confident regions (hap.py specific)
            return True
        if not happy_vcf and var.filter.get('PASS'):
            # Don't check non-PASS variants
            return True
        return False

    for var in f_vcf.fetch():
        if select_variant(var):
            annotated_v = allele_match.fetch_nearby_cohort(
                var=var, f_query_vcf=f_query_vcf,
                f_fasta=f_fasta, 
                update_info=update_info,
                query_info=None,
                debug=debug)
            if annotated_v:
                f_out.write(annotated_v)

if __name__ == '__main__':
    # We use python3.6 because we require dict is ordered.
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()

    match_vcf(
        fn_vcf=args.vcf,
        fn_query_vcf=args.query_vcf,
        fn_fasta=args.ref,
        fn_out=args.out,
        happy_vcf=args.happy,
        debug=args.debug
    )
