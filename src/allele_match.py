import argparse
import pysam
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='path to target happy.out.vcf'
    )
    parser.add_argument(
        '-p', '--panel',
        help='path to the reference panel VCF (TBI or CSI indexes are required)'
    )
    parser.add_argument(
        '-r', '--ref',
        help='path to the reference FASTA (FAI index is required)'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='path to output VCF. Set to "-" to print to stdout. ["-"]'
    )
    args = parser.parse_args()
    return args


def compare_haplotypes(
    var: pysam.VariantRecord, cohort_vars: list, ref: str
) -> pysam.VariantRecord:
    # print(var)
    # for c in cohort_vars:
    #     print(c)
    #     print(c.contig)
    #     print(c.start)
    #     print(c.stop)
    #     print(c.alleles)
    # print('ref', ref)
    start = min(var.start, min([v.start for v in cohort_vars]))
    # print(var.contig)
    # print(var.start)
    # print(var.stop)
    # print(var.alleles)

    # dict_alt_af:
    #   - key: local haplotype
    #   - value: allele frequency
    dict_alt_af = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_alt_af[var_seq] = 0
        # print(dict_alt_af)
    for c_var in cohort_vars:
        # loop through matched cohort variants
        for i, alt in enumerate(c_var.alts):
            # loop through each allele (index=`i`) in a cohort variant
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            # print(c_var_seq, c_var.info['AF'][i])
            if c_var_seq in dict_alt_af:
                try:
                    dict_alt_af[c_var_seq] = c_var.info['AF'][i]
                except:
                    raise ValueError('Error: "AF" field is not provided in a cohort variant')
    var.info.__setitem__('AF', tuple(dict_alt_af.values()))
    
    return var


def match_allele(
    var: pysam.VariantRecord,
    f_panel: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    f_out: pysam.VariantFile
) -> None:
    # var.start: 0-based; var.pos: 1-based
    # Pysam uses 0-based
    # var_region = (var.contig, var.start, var.start + max(var.alleles))
    # Fetch cohort variants
    var_maxstop = max([var.start + len(a) for a in var.alleles])
    cohort_vars = list(f_panel.fetch(
        var.contig, var.start, var_maxstop))

    if len(cohort_vars) == 0:
        # If cannot find matched cohorts, set AF to 0
        var.info.__setitem__('AF', 0)
    else:
        cohort_start = min(var.start, min([v.start for v in cohort_vars]))
        cohort_maxstop = var_maxstop
        for v in cohort_vars:
            cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) for a in v.alleles]))

        if not all([v.contig == var.contig for v in cohort_vars]):
            # All variants should have the same contig
            raise ValueError(
                "Fetched variants have disconcordant contigs: ",
                [v.contig for v in cohort_vars])
        # Fetch reference sequence
        try:
            ref_seq = f_fasta.fetch(reference=var.contig, start=cohort_start, end=cohort_maxstop)
        except:
            var.info.__setitem__('AF', 0)
            print('Warning: encounter the edge of a contig. Set "AF"=0', file=sys.stderr)
            # raise ValueError("Errors during fetching allele matching sequence in the ref FASTA")
        f_out.write(compare_haplotypes(var, cohort_vars, ref_seq))


def annotate_vcf(
    fn_vcf: str, fn_panel_vcf: str, fn_fasta: str, fn_out: str
) -> None:
    try:
        f_vcf = pysam.VariantFile(fn_vcf)
        f_vcf.header.add_meta('INFO', items=[('ID','AF'), ('Number','A'), ('Type','Float'), ('Description','Population allele frequency')])
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
    
    for var in f_vcf.fetch():
        if len(var.filter.keys()) != 1:
            print('Error: more than one filters for a variant. Exit.', file=sys.stderr)
            print(var)
            exit(1)
        # Only check variants in confident regions (hap.py specific)
        # if var.info.get('Regions'):
        elif var.filter.keys()[0] == 'PASS':
            # Only take 'PASS' variants
            match_allele(var, f_panel, f_fasta, f_out)


if __name__ == '__main__':
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()

    annotate_vcf(
        fn_vcf=args.vcf,
        fn_panel_vcf=args.panel,
        fn_fasta=args.ref,
        fn_out=args.out
    )