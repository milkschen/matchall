'''
An haplotype-based matching algorithm to match alleles even when they are in different forms.
'''
import pysam
import sys


def update_info_empty(var: pysam.VariantRecord, info: dict):
    info_num = info['Number']
    if info_num == 'A':
        var.info.__setitem__(info['ID'], tuple([0 for _ in var.alts]))
    # Legal format, but have not supported
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
    update_info: dict,
    query_info: dict=None,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' 
    Match a variant with nearby cohorts using local haplotypes.

    Inputs:
        - var: Target variant.
        - cohort_vars: list of pysam.VariantRecord. Fetched nearby cohorts.
        - ref: Local REF haplotype
        - update_info: VCF INFO field to update. 'ID', 'Number', 'Type' are required.
            E.g.
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'}
        - query_info: VCF INFO field to query. If not set, use `update_info`.

    Returns:
        - var: Target variant with annotation.
    
    Raises:
        - ValueError: 
            * If query_info['ID'] is not found in a queried variant
            * If pysam.VariantRecord.__setitem__() fails
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
    
    # if not query_info:
    #     query_info = update_info
    update_info_num = update_info['Number']
    if query_info:
        query_info_num = query_info['Number']
    # dict_hap_info:
    #   - key: local haplotype
    #   - value: queried info value
    dict_hap_info = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_hap_info[var_seq] = 0

    # Loop through matched cohort variants
    for c_var in cohort_vars:
        # print(c_var.info[info['ID']])
        # Loop through each allele (index=`i`) in a cohort variant
        for i, alt in enumerate(c_var.alts):
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            if c_var_seq in dict_hap_info:
                try:
                    if not query_info:
                        dict_hap_info[c_var_seq] = 1
                    elif query_info_num in ['0', '1']:
                        dict_hap_info[c_var_seq] = c_var.info[query_info['ID']]
                    else:
                        dict_hap_info[c_var_seq] = c_var.info[query_info['ID']][i]
                except:
                    raise ValueError(f'Error: INFO."{query_info["ID"]}" is in the fetched cohort variant')
    
    try:
        if update_info_num in ['0', '1']:
            var.info.__setitem__(update_info['ID'], tuple(dict_hap_info.values())[0])
        # Rare weird cases where both alts are the same
        elif len(dict_hap_info.keys()) != len(var.alts):
            var.info.__setitem__(update_info['ID'], tuple([list(dict_hap_info.values())[0] for _ in var.alts]))
        else:
            var.info.__setitem__(update_info['ID'], tuple(dict_hap_info.values()))
    except:
        raise ValueError(f'Error: VariantRecord.__setitem__ failed')
    
    return var


def fetch_nearby_cohort(
    var: pysam.VariantRecord,
    f_query_vcf: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    update_info: dict,
    query_info: dict=None,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Fetch nearby cohorts and local REF haplotype for a variant.

    Inputs:
        - var: Target variant.
        - f_query_vcf: Queried VCF file (e.g. a cohort VCF).
        - f_fasta: REF FASTA file.
        - update_info: VCF INFO field to update. 'ID', 'Number', 'Type' are required.
            E.g.
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'}
        - query_info: VCF INFO field to query. If not set, use `update_info`.

    Raises:
        - ValueError: If fetched variants don't share the same contig.
    '''
    # var.start: 0-based; var.pos: 1-based
    # Pysam uses 0-based
    # var_region = (var.contig, var.start, var.start + max(var.alleles))
    # Fetch cohort variants
    var_maxstop = max([var.start + len(a) for a in var.alleles])
    cohort_vars = list(f_query_vcf.fetch(
        var.contig, var.start, var_maxstop))

    # If cannot find matched cohorts, set AF to 0
    if len(cohort_vars) == 0:
        update_info_empty(var, update_info)
        return var
    
    cohort_start = min(var.start, min([v.start for v in cohort_vars]))
    cohort_maxstop = var_maxstop
    for v in cohort_vars:
        cohort_maxstop = max(cohort_maxstop,
                                max([v.start + len(a) for a in v.alleles]))

    # All variants should have the same contig
    if not all([v.contig == var.contig for v in cohort_vars]):
        raise ValueError(
            "Fetched variants have disconcordant contigs: ",
            [v.contig for v in cohort_vars])
    
    # Fetch reference sequence
    try:
        ref_seq = f_fasta.fetch(
            reference=var.contig, start=cohort_start, end=cohort_maxstop)
    except:
        update_info_empty(var, update_info)
        print(f'Warning: encounter the edge of a contig. Set "{update_info["ID"]}" as the init value.', file=sys.stderr)
    
    try:
        return match_allele(
            var=var, cohort_vars=cohort_vars, 
            ref=ref_seq, 
            update_info=update_info,
            query_info=query_info,
            debug=debug)
    except Exception as e:
        print(f'Error: {e}')
        return None
