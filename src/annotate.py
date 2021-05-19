'''
Annotate a VCF using info from another VCF.

An haplotype-based matching algorithm is used to match alleles represented in different forms.

Example:
python annotate.py -v <target.vcf.gz> -p <panel.vcf.gz> -r <ref.fa> -o <target.annotated.vcf.gz>
'''

import allele_match
import argparse
import pysam
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [required]'
    )
    parser.add_argument(
        '-q', '--query-vcf', required=True,
        help='Path to query VCF (TBI or CSI indexes are required). This can be a reference panel. [required]'
    )
    parser.add_argument(
        '-r', '--ref',
        help='Path to the reference FASTA (FAI index is required). [required]'
    )
    parser.add_argument(
        '--info', default='AF',
        help='INFO tag to query. ["AF"]'
    )
    # parser.add_argument(
    #     '--info-description', default='Population allele frequency',
    #     help='Description of the INFO tag to query. ["Population allele frequency"]'
    # )
    # parser.add_argument(
    #     '--allele-frequency-cutoff', type=float, default=0,
    #     help= \
    #         'allele frequency cutoff (set to a non-zero value to activate).\n' +
    #         'Split output VCFs will be write to ' +
    #         '`allele-frequency-prefix`-af_gt_`allele-frequency-cutoff`.vcf' +
    #         'and `allele-frequency-prefix`-af_leq_`allele-frequency-cutoff`.vcf\n' +
    #         '`allele-frequency-prefix` must be set, too [0]'
    # )
    # parser.add_argument(
    #     '--allele-frequency-prefix', default=None,
    #     help= \
    #         'Prefix of output files in the `allele-frequency-cutoff` mode.\n' +
    #         'See `allele-frequency-cutoff` for more explanation. [None]'
    # )
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
    fn_vcf: str,
    fn_query_vcf: str,
    fn_fasta: str,
    fn_out: str,
    info_tag: str,
    happy_vcf: bool,
    debug: bool=False
    # af_cutoff: float=0, af_prefix: str=None
    ) -> None:
    # Open panel VCF. Searching for `info_tag` in the header
    info = None
    try:
        f_query_vcf = pysam.VariantFile(fn_query_vcf)
        for r in f_query_vcf.header.records:
            if r.type == 'INFO':
                if r.get('ID') == info_tag:
                    info = dict(r)
                    break    
    except:
        raise ValueError(f'Error: Cannot open "{fn_query_vcf}"')

    # If info_tag cannot be found, exit program
    if info is None:
        print(f'Error: info-tag "{info_tag}" cannot be found in f{fn_query_vcf}. Exit.')
        exit(1)
    info.pop('IDX', None)
    info['Description'] = info['Description'].replace('"', '')
    
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        f_vcf.header.add_meta('INFO', items = info.items())
    except:
        raise ValueError(f'Error: Cannot open "{fn_vcf}"')
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
                update_info=info,
                query_info=info,
                debug=debug)
            if annotated_v:
                f_out.write(annotated_v)


if __name__ == '__main__':
    # We use python3.6 because we require dict is ordered.
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()

    annotate_vcf(
        fn_vcf=args.vcf,
        fn_query_vcf=args.query_vcf,
        fn_fasta=args.ref,
        fn_out=args.out,
        info_tag=args.info,
        happy_vcf=args.happy,
        debug=args.debug
        # af_cutoff=args.allele_frequency_cutoff,
        # af_prefix=args.allele_frequency_prefix
    )
