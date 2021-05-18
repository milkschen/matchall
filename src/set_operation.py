'''
Perform set operations on two VCFs.
E.g. Finding the intersection/union set.
'''

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
        '-o', '--out', required=True,
        help='Prefix to output VCF. [required]'
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


def match_allele(
    var: pysam.VariantRecord, cohort_vars: list, 
    ref: str, 
    # info_field: str, 
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Match a variant with nearby cohorts using local haplotypes.

    Inputs:
        - var: Target variant.
        - cohort_vars: list of pysam.VariantRecord. Fetched nearby cohorts.
        - ref: Local REF haplotype
        # - info_field: the INFO field interested in querying
    Returns:
        - var: Target variant with annotation.
    Raises:
        - ValueError: If `info_field` is not found in cohort variants.
    '''
    if debug:
        print('var')
        print(var)
        print(var.start, var.stop, var.alleles, var.alts)
        print('cohorts')
        for c in cohort_vars:
            print(c)
            print(c.start, c.stop, c.alleles)
        print('ref')
        print(ref)
    start = min(var.start, min([v.start for v in cohort_vars]))

    # dict_alt:
    #   - key: local haplotype
    #   - value: info_field
    dict_alt = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_alt[var_seq] = 0
    for c_var in cohort_vars:
        # loop through matched cohort variants
        for i, alt in enumerate(c_var.alts):
            # loop through each allele (index=`i`) in a cohort variant
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            if c_var_seq in dict_alt:
                try:
                    dict_alt[c_var_seq] = 1
                    # dict_alt[c_var_seq] = c_var.info[info_field][i]
                except:
                    raise ValueError()
                    # raise ValueError(
                    #     f'Error: "{info_field}" field is not provided in a cohort variant')

    if len(dict_alt.keys()) != len(var.alts):
        # Rare weird cases where both alts are the same
        # print('no_eq')
        # print(tuple([list(dict_alt.values())[0] for i in var.alts]))
        var.info.__setitem__('MATCH', tuple([list(dict_alt.values())[0] for i in var.alts]))
    else:
        # print('eq')
        # print(var)
        # print(tuple(dict_alt.values()))
        var.info.__setitem__('MATCH', tuple(dict_alt.values()))
        # print(var)
    
    return var


def fetch_nearby_cohort(
    var: pysam.VariantRecord,
    f_panel: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    f_out: pysam.VariantFile,
    info_field: str, 
    debug: bool=False
    ) -> None:
    ''' Fetch nearby cohorts and local REF haplotype for a variant.

    Inputs:
        - var: Target variant.
        - f_panel: Cohort VCF file.
        - f_fasta: REF FASTA file.
        - f_out: Output VCF file.
    Raises:
        - ValueError: If fetched variants don't share the same contig.
    '''
    # var.start: 0-based; var.pos: 1-based
    # Pysam uses 0-based
    # var_region = (var.contig, var.start, var.start + max(var.alleles))
    # Fetch cohort variants
    var_maxstop = max([var.start + len(a) for a in var.alleles])
    cohort_vars = list(f_panel.fetch(
        var.contig, var.start, var_maxstop))

    if len(cohort_vars) == 0:
        # If cannot find matched cohorts, set `info_field` to 0
        var.info.__setitem__(info_field, tuple([0 for i in var.alts]))
        f_out.write(var)
    else:
        cohort_start = min(var.start, min([v.start for v in cohort_vars]))
        cohort_maxstop = var_maxstop
        for v in cohort_vars:
            cohort_maxstop = max(cohort_maxstop,
                                 max([v.start + len(a) for a in v.alleles]))

        if not all([v.contig == var.contig for v in cohort_vars]):
            # All variants should have the same contig
            raise ValueError(
                "Fetched variants have disconcordant contigs: ",
                [v.contig for v in cohort_vars])
        # Fetch reference sequence
        try:
            ref_seq = f_fasta.fetch(
                reference=var.contig, start=cohort_start, end=cohort_maxstop)
        except:
            var.info.__setitem__(info_field, tuple([0 for i in var.alts]))
            print('Warning: encounter the edge of a contig. Set "AF"=0', file=sys.stderr)
        
        try:
            f_out.write(match_allele(var, cohort_vars, ref_seq, debug))
        except:
            print('Warning: unexpected error when fetching nearby cohort variants')


def annotate_vcf(
    fn_vcf: str, fn_query_vcf: str, fn_fasta: str,
    fn_out: str, happy_vcf: bool, debug: bool=False
    # af_cutoff: float=0, af_prefix: str=None
    ) -> None:
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        f_vcf.header.add_meta(
            'INFO',
            items=[('ID','AF'), ('Number','A'),
                   ('Type','Float'), ('Description','Population allele frequency')])
        f_vcf.header.add_meta(
            'INFO',
            items=[('ID','MATCH'), ('Number','A'),
                   ('Type','Integer'), ('Description','If genotype is matched with a query')])
    except:
        raise ValueError(f'Error: Cannot open "{fn_vcf}"')
    try:
        f_panel = pysam.VariantFile(fn_query_vcf)
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

    for var in f_vcf.fetch():
        if happy_vcf:
            # Only check variants in confident regions (hap.py specific)
            if var.info.get('Regions'):
                fetch_nearby_cohort(
                    var=var, f_panel=f_panel, f_fasta=f_fasta, 
                    f_out=f_out, info_field='MATCH', debug=debug)
        else:
            if var.filter.get('PASS'):
                # Only take 'PASS' variants
                fetch_nearby_cohort(
                    var=var, f_panel=f_panel, f_fasta=f_fasta, 
                    f_out=f_out, info_field='MATCH', debug=debug)

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
        happy_vcf=args.happy,
        debug=args.debug
        # af_cutoff=args.allele_frequency_cutoff,
        # af_prefix=args.allele_frequency_prefix
    )
