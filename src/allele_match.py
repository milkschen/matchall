'''
An haplotype-based matching algorithm to match alleles even when they are in different forms.
'''
import pysam
import sys


def match_allele(
    var: pysam.VariantRecord, cohort_vars: list, ref: str, debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Match a variant with nearby cohorts using local haplotypes.

    Inputs:
        - var: Target variant.
        - cohort_vars: list of pysam.VariantRecord. Fetched nearby cohorts.
        - ref: Local REF haplotype
    Returns:
        - var: Target variant with annotation.
    Raises:
        - ValueError: If "AF" is not found in cohort variants.
    '''
    if debug:
        print(f'var = {var}'.rstrip())
        print(f'    start={var.start}')
        print(f'    stop={var.stop}')
        print(f'    alleles={var.alleles}')
        print(f'    alts={var.alts}')
        print('cohorts')
        for c in cohort_vars:
            print(f'{c}'.rstrip())
            print('   ', c.start, c.stop, c.alleles)
        print('ref =', ref)
    start = min(var.start, min([v.start for v in cohort_vars]))

    # dict_alt_af:
    #   - key: local haplotype
    #   - value: allele frequency
    dict_alt_af = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_alt_af[var_seq] = 0
    for c_var in cohort_vars:
        # loop through matched cohort variants
        for i, alt in enumerate(c_var.alts):
            # loop through each allele (index=`i`) in a cohort variant
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            if c_var_seq in dict_alt_af:
                try:
                    dict_alt_af[c_var_seq] = c_var.info['AF'][i]
                except:
                    raise ValueError('Error: "AF" field is not provided in a cohort variant')

    if len(dict_alt_af.keys()) != len(var.alts):
        # Rare weird cases where both alts are the same
        var.info.__setitem__('AF', tuple([list(dict_alt_af.values())[0] for i in var.alts]))
    else:
        var.info.__setitem__('AF', tuple(dict_alt_af.values()))
    
    return var


def fetch_nearby_cohort(
    var: pysam.VariantRecord,
    f_panel: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Fetch nearby cohorts and local REF haplotype for a variant.

    Inputs:
        - var: Target variant.
        - f_panel: Cohort VCF file.
        - f_fasta: REF FASTA file.
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

    # If cannot find matched cohorts, set AF to 0
    if len(cohort_vars) == 0:
        var.info.__setitem__('AF', tuple([0 for i in var.alts]))
        return var
    
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
        var.info.__setitem__('AF', tuple([0 for i in var.alts]))
        print('Warning: encounter the edge of a contig. Set "AF"=0', file=sys.stderr)
        # raise ValueError("Errors during fetching allele matching sequence in the ref FASTA")
    
    try:
        return match_allele(var, cohort_vars, ref_seq, debug)
    except:
        print('Warning: unexpected error at allele_match.py:match_allele()')
        return None
