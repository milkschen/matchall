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


def match_vcf(
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
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        # f_vcf.header.info.get('INFO', update_info.items())
        f_vcf.header.add_meta('INFO', items=update_info.items())
    except Exception as e:
        raise(e)
    # except:
    #     raise ValueError(f'Error: Cannot open "{fn_vcf}"')

    try:
        f_query_vcf = pysam.VariantFile(fn_query_vcf)
        if do_isec or do_private:
            f_query_vcf.header.add_meta('INFO', items=update_info.items())
    except Exception as e:
        raise(e)
        # raise ValueError(f'Error: Cannot open "{fn_query_vcf}"')

    try:
        f_fasta = pysam.FastaFile(fn_fasta)
    except:
        raise ValueError(f'Error: Cannot open "{fn_fasta}"')

    if do_annotate: 
        try:
            f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
        except:
            raise ValueError(f'Error: Cannot create "{fn_out}"')
    
    if do_isec: 
        try:
            fn_isec0 = out_prefix + '.isec0.vcf.gz'
            f_isec0 = pysam.VariantFile(fn_isec0, 'w', header=f_vcf.header)
        except:
            raise ValueError(f'Error: Cannot create "{fn_isec0}"')
        try:
            fn_isec1 = out_prefix + '.isec1.vcf.gz'
            f_isec1 = pysam.VariantFile(fn_isec1, 'w', header=f_vcf.header)
        except:
            raise ValueError(f'Error: Cannot create "{fn_isec1}"')

    if do_private: 
        try:
            fn_private0 = out_prefix + '.private0.vcf.gz'
            f_private0 = pysam.VariantFile(fn_private0, 'w', header=f_vcf.header)
        except:
            raise ValueError(f'Error: Cannot create "{fn_private0}"')
        try:
            fn_private1 = out_prefix + '.private1.vcf.gz'
            f_private1 = pysam.VariantFile(fn_private1, 'w', header=f_vcf.header)
        except:
            raise ValueError(f'Error: Cannot create "{fn_private1}"')

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
            try:
                if do_annotate and annotated_v:
                    f_out.write(annotated_v)
                
                if do_isec and annotated_v:
                    # if any(annotated_v.info['MATCH']):
                    if annotated_v.info['MATCH']:
                        f_isec0.write(annotated_v)
                
                if do_private and annotated_v:
                    # if not all(annotated_v.info['MATCH']):
                    if not annotated_v.info['MATCH']:
                        f_private0.write(annotated_v)
            except Exception as e:
                print(e)
    
    # Second loop - using f_query_vcf as main and f_vcf as query
    if do_isec or do_private:
        for var in f_query_vcf.fetch():
            if select_variant(var):
                annotated_v = allele_match.fetch_nearby_cohort(
                    var=var, f_query_vcf=f_vcf,
                    f_fasta=f_fasta, 
                    update_info=update_info,
                    query_info=None,
                    debug=debug)
                
                try:
                    if do_isec and annotated_v:
                        # if any(annotated_v.info['MATCH']):
                        if annotated_v.info['MATCH']:
                            f_isec1.write(annotated_v)
                    
                    if do_private and annotated_v:
                        # if not all(annotated_v.info['MATCH']):
                        if not annotated_v.info['MATCH']:
                            f_private1.write(annotated_v)
                except Exception as e:
                    print(e)


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
        out_prefix=args.out_prefix,
        mode=args.mode,
        happy_vcf=args.happy,
        debug=args.debug
    )
