'''
An haplotype-based matching algorithm to match alleles even when they are in different forms.
'''
import pysam
import sys


def update_info_empty(var: pysam.VariantRecord, info: dict):
    info_num = info['Number']
    if info_num == 'A':
        var.info.__setitem__(info['ID'], tuple([0 for _ in var.alts]))
    # elif info_num == 'R':
    #     var.info.__setitem__(info['ID'], tuple([0 for _ in var.alleles]))
    elif info_num == '0':
        var.info.__setitem__(info['ID'], False)
    elif info_num == '1':
        var.info.__setitem__(info['ID'], 0)
    elif info_num.isnumeric():
        var.info.__setitem__(info['ID'], tuple([0 for _ in int(info_num)]))
    else:
        raise(ValueError, f'Unsupported format: "{info_num}"')


def match_allele(
    var: pysam.VariantRecord,
    cohort_vars: list,
    ref: str,
    info: dict,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' 
    Match a variant with nearby cohorts using local haplotypes.

    Inputs:
        - var: Target variant.
        - cohort_vars: list of pysam.VariantRecord. Fetched nearby cohorts.
        - ref: Local REF haplotype
        - info: VCF INFO field. 'ID', 'Number', 'Type' are required.
            E.g.
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'}

    Returns:
        - var: Target variant with annotation.
    
    Raises:
        - ValueError: 
            * If `info_tag` is not found in cohort variants
            * If __setitem__ fails
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
    
    info_num = info['Number']
    # dict_alt_info:
    #   - key: local haplotype
    #   - value: queried `info_tag` value
    dict_alt_info = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_alt_info[var_seq] = 0

    # Loop through matched cohort variants
    for c_var in cohort_vars:
        # print(c_var.info[info['ID']])
        # Loop through each allele (index=`i`) in a cohort variant
        for i, alt in enumerate(c_var.alts):
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            if c_var_seq in dict_alt_info:
                try:
                    if info_num in ['0', '1']:
                        dict_alt_info[c_var_seq] = c_var.info[info['ID']]
                    else:
                        dict_alt_info[c_var_seq] = c_var.info[info['ID']][i]
                except:
                    raise ValueError(f'Error: INFO."{info["ID"]}" is not provided in the fetched cohort variant')
    
    try:
        if info_num in ['0', '1']:
            var.info.__setitem__(info['ID'], tuple(dict_alt_info.values())[0])
        elif len(dict_alt_info.keys()) != len(var.alts):
            # Rare weird cases where both alts are the same
            var.info.__setitem__(info['ID'], tuple([list(dict_alt_info.values())[0] for _ in var.alts]))
        else:
            var.info.__setitem__(info['ID'], tuple(dict_alt_info.values()))
    except:
        raise ValueError(f'Error: VariantRecord.__setitem__ failed')
    
    return var


def fetch_nearby_cohort(
    var: pysam.VariantRecord,
    f_panel: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    info: dict,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Fetch nearby cohorts and local REF haplotype for a variant.

    Inputs:
        - var: Target variant.
        - f_panel: Cohort VCF file.
        - f_fasta: REF FASTA file.
        - info: VCF INFO field. 'ID', 'Number', 'Type' are required.
            E.g.
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'}

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
        update_info_empty(var, info)
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
        update_info_empty(var, info)
        print(f'Warning: encounter the edge of a contig. Set "{info["ID"]}" as the init value.', file=sys.stderr)
        # raise ValueError("Errors during fetching allele matching sequence in the ref FASTA")
    
    try:
        return match_allele(
            var=var, cohort_vars=cohort_vars, 
            ref=ref_seq, 
            info=info,
            debug=debug)
    except:
        print('Warning: unexpected error at allele_match.py:match_allele()')
        return None
