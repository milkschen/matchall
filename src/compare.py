'''
Compare two VCFs using the matchall algorithm.
This script supports the below operations:
    - [annotate] annotate matching status of variants using a "MATCH" tag in the INFO field
    - [isec] find the intersected variants with respect to each input VCF
    - [private] find variant sets private to each input VCF

Usage:
    python src/compare.py -v A.vcf.gz -q B.vcf.gz -op A_0-B_1 -m annotate,private,isec \
        -o A_0-B_1.vcf.gz -r genome.fa
'''
import matchall
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
        '-op', '--out-prefix', default=None,
        help='Prefix to output VCF for "isec" and "private" modes. [None]'
    )
    parser.add_argument(
        '-m', '--mode', default='annotate',
        help='Mode: ["annotate", "isec", "private"]. ' + 
             'Multiple modes can be toggled at once, e.g. "annotate,private,isec". ["annotate"]'
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


def write_to_isec_and_private(var, do_isec, do_private, f_isec, f_private) -> None:
    ''' Write a variant to private/isec VCFs.

    Criteria:
        - if all matched -> isec
        - if none matched -> private
        - mixed: 
            * matched alleles -> isec (unmatched alleles masked using "*")
            * unmathced alleles -> private (matched alleles masked using "*")

    '''
    if (not do_isec) and (not do_private):
        return
    
    match_status = var.info['MATCH']
    if all(match_status) == 1:
        if do_isec:
            f_isec.write(var)
    elif not any(match_status) == 1:
        if do_private:
            f_private.write(var)
    else: # compound
        for i, m in enumerate(match_status):
            if m and do_isec:
                alleles = list(var.alleles)
                alleles[i + 1] = '*'
                var.alleles = alleles
                alts = list(var.alts)
                alts[i] = '*'
                var.alts = alts
                if not all([a == '*' for a in alts]):
                    f_isec.write(var)
            elif (not m) and do_private:
                alleles = list(var.alleles)
                alleles[i + 1] = '*'
                var.alleles = alleles
                alts = list(var.alts)
                alts[i] = '*'
                var.alts = alts
                if not all([a == '*' for a in alts]):
                    f_private.write(var)


def compare_vcf(
    fn_vcf: str,
    fn_query_vcf: str,
    fn_fasta: str,
    fn_out: str,
    mode: list,
    out_prefix: str=None,
    happy_vcf: bool=False,
    debug: bool=False
    ) -> None:
    mode = mode.split(',')
    do_annotate, do_isec, do_private = False, False, False
    for m in mode:
        if m not in ['annotate', 'isec', 'private']:
            raise(ValueError, f'Error: unknown mode {m}')
        if m == 'annotate':
            do_annotate = True
        elif m == 'isec':
            do_isec = True
        elif m == 'private':
            do_private = True
    if do_isec or do_private:
        if out_prefix is None:
            raise(
                ValueError,
                f'Error: `output_prefix` must be set under "isec" and "private" modes.')

    update_info = {'ID': 'MATCH', 'Number': 'A', 'Type': 'Integer', 'Description': 'If genotype is matched with a query'}
    # update_info = {'ID': 'MATCH', 'Number': '1', 'Type': 'Integer', 'Description': 'If genotype is matched with a query'}

    f_vcf = pysam.VariantFile(fn_vcf)
    f_vcf.header.add_meta('INFO', items=update_info.items())
    f_query_vcf = pysam.VariantFile(fn_query_vcf)
    if do_isec or do_private:
        f_query_vcf.header.add_meta('INFO', items=update_info.items())
    f_fasta = pysam.FastaFile(fn_fasta)
    
    if do_annotate: 
        f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
    
    if do_isec: 
        fn_isec0 = out_prefix + '.isec0.vcf.gz'
        f_isec0 = pysam.VariantFile(fn_isec0, 'w', header=f_vcf.header)
        fn_isec1 = out_prefix + '.isec1.vcf.gz'
        f_isec1 = pysam.VariantFile(fn_isec1, 'w', header=f_query_vcf.header)

    if do_private: 
        fn_private0 = out_prefix + '.private0.vcf.gz'
        f_private0 = pysam.VariantFile(fn_private0, 'w', header=f_vcf.header)
        fn_private1 = out_prefix + '.private1.vcf.gz'
        f_private1 = pysam.VariantFile(fn_private1, 'w', header=f_query_vcf.header)

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
            try:
                annotated_v = matchall.fetch_nearby_cohort(
                    var=var, f_query_vcf=f_query_vcf,
                    f_fasta=f_fasta, 
                    update_info=update_info,
                    query_info=None,
                    debug=debug)
                if do_annotate and annotated_v:
                    f_out.write(annotated_v)
                
                write_to_isec_and_private(
                    annotated_v, do_isec, do_private, f_isec0, f_private0)

            except Exception as e:
                print(f'Warning: encounter the below exception when querying {fn_vcf} agains {fn_query_vcf}')
                print(e)
                print(var)
    
    # Second loop - using f_query_vcf as main and f_vcf as query
    if do_isec or do_private:
        for var in f_query_vcf.fetch():
            if select_variant(var):
                try:
                    annotated_v = matchall.fetch_nearby_cohort(
                        var=var, f_query_vcf=f_vcf,
                        f_fasta=f_fasta, 
                        update_info=update_info,
                        query_info=None,
                        debug=debug)
                
                    write_to_isec_and_private(
                        annotated_v, do_isec, do_private, f_isec1, f_private1)
                    
                except Exception as e:
                    print(f'Warning: encounter the below exception when querying {fn_query_vcf} agains {fn_vcf}')
                    print(e)
                    print(var)


if __name__ == '__main__':
    # We use python3.6 because we require dict is ordered.
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()

    compare_vcf(
        fn_vcf=args.vcf,
        fn_query_vcf=args.query_vcf,
        fn_fasta=args.ref,
        fn_out=args.out,
        out_prefix=args.out_prefix,
        mode=args.mode,
        happy_vcf=args.happy,
        debug=args.debug
    )
