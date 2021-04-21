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
        help='path to output VCF. Set to "-" to print to stdout. ["-"]'
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


def compare_haplotypes(
    var: pysam.VariantRecord, cohort_vars: list, ref: str
) -> pysam.VariantRecord:
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
        var.info.__setitem__('AF', tuple([0 for i in var.alts]))
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
            var.info.__setitem__('AF', tuple([0 for i in var.alts]))
            print('Warning: encounter the edge of a contig. Set "AF"=0', file=sys.stderr)
            # raise ValueError("Errors during fetching allele matching sequence in the ref FASTA")
        f_out.write(compare_haplotypes(var, cohort_vars, ref_seq))


def annotate_vcf(
    fn_vcf: str, fn_panel_vcf: str, fn_fasta: str, fn_out: str, happy_vcf: bool,
    # af_cutoff: float=0, af_prefix: str=None
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

    # if af_cutoff > 0 and af_prefix == None:
    #     raise ValueError(f'Error: `allele-frequency-prefix` needs to be set when `allele-frequency-cutoff` > 0')
    # elif af_cutoff > 0:
    #     f_out_high = pysam.VariantFile(af_prefix+f'-af_gt_{af_cutoff}.vcf', 'w', header=f_vcf.header)
    #     f_out_low = pysam.VariantFile(af_prefix+f'-af_leq_{af_cutoff}.vcf', 'w', header=f_vcf.header)
    
    for var in f_vcf.fetch():
        # Only check variants in confident regions (hap.py specific)
        if happy_vcf:
            if var.info.get('Regions'):
                match_allele(var, f_panel, f_fasta, f_out)
        else:
            # if len(var.filter.keys()) != 1:
            #     print('Warning: more than one filters for a variant. Exit.', file=sys.stderr)
            #     print(var)
            #     continue
            if var.filter.get('PASS'):
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
        fn_out=args.out,
        happy_vcf=args.happy
        # af_cutoff=args.allele_frequency_cutoff,
        # af_prefix=args.allele_frequency_prefix
    )
