import argparse
import os
import sys
import pandas as pd
import pysam


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
        '-f', '--fasta',
        help='path to the reference FASTA (FAI index is required)'
    )
    parser.add_argument(
        '-o', '--out',
        help='path to output TSV'
    )
    parser.add_argument(
        '-f', '--force-run', action='store_true',
        help='Set to force running [False]'
    )
    args = parser.parse_args()
    return args


def update_happy_outcome(var, outcome):
    bvt = []
    bd = []
    # happy samples: 'TRUTH' and 'QUERY'
    for sample in var.samples.items():
        for v_format in sample[1].items():
            if v_format[0] == 'BD':
                bd.append(v_format[1])
            elif v_format[0] == 'BVT':
                bvt.append(v_format[1])
    # TRUTH
    if bd[0] == 'TP':
        outcome[bvt[0]]['TRUTH.TP'] += 1
        outcome[bvt[0]]['TRUTH.TOTAL'] += 1
    elif bd[0] == 'FN':
        outcome[bvt[0]]['TRUTH.FN'] += 1
        outcome[bvt[0]]['TRUTH.TOTAL'] += 1
    # QUERY
    if bd[1] == 'FP':
        outcome[bvt[1]]['QUERY.FP'] += 1


def write_to_tsv(outcome, out_fn):
    print(outcome, file=sys.stderr)
    df = pd.DataFrame.from_dict(outcome, orient='index')
    df['METRIC.Recall'] = df['TRUTH.TP'] / (df['TRUTH.TP'] + df['TRUTH.FN'])
    df['METRIC.Precision'] = df['TRUTH.TP'] / (df['TRUTH.TP'] + df['QUERY.FP'])
    df['METRIC.F1_Score'] = 2 / ((1/df['METRIC.Recall']) + (1/df['METRIC.Precision']))
    df.to_csv(out_fn, sep='\t', float_format='%.6f')
    print(df, file=sys.stderr)


def summarize_accuracy(fn_happy_vcf, fn_summary):
    ''' Summarize the accuracy of a hap.py annotated VCF file. '''
    cycle = 0
    fields = {
        'TRUTH.TOTAL': 0,
        'TRUTH.TP': 0,
        'TRUTH.FN': 0,
        'QUERY.FP': 0,
    }
    outcome = {'SNP': fields, 'INDEL': fields.copy()}

    f_happy_vcf = pysam.VariantFile(fn_happy_vcf)
    for var in f_happy_vcf.fetch():
        # Only check variants in confident regions
        if var.info.get('Regions'):
            update_happy_outcome(var, outcome)
            
            if cycle > 0 and cycle % 1000000 == 0:
                print(f'Process {cycle} records', file=std.stderr)
            cycle += 1

    write_to_tsv(outcome, fn_summary)



def match_allele(var, f_panel, f_fasta):
    # var.start: 0-based; var.pos: 1-based
    # Pysam uses 0-based
    # var_region = (var.contig, var.start, var.start + max(var.alleles))
    # Fetch cohort variants
    cohort_vars = f_panel.fetch(var.contig, var.start, var.start + max(var.alleles))
    cohort_start = min([v.start for v in cohort_vars])
    cohort_stop = min([v.stop for v in cohort_vars])
    if not all([v.contig == var.contig for v in cohort_vars]):
        raise ValueError(
            "Fetched variants have disconcordant contigs: ",
            [v.contig for v in cohort_vars])
    # Fetch reference sequence
    try:
        query_sequence = f_fasta.fetch(var.contig, cohort_start, cohort_stop)
    except:
        raise ValueError("Errors during fetching allele matching sequence in the ref FASTA")
    

    pass


def annotate_vcf(fn_vcf, fn_panel_vcf, fn_fasta, fn_summary):
    f_vcf = pysam.VariantFile(fn_vcf)
    f_panel = pysam.VariantFile(fn_panel_vcf)
    f_fasta = pysam.FastaFile
    # for filt in f_vcf.header.filters:
    #     print(filt)
    # cnt = 0
    for var in f_vcf.fetch():
        # hap.py specific
        # # Only check variants in confident regions
        # if var.info.get('Regions'):
        #     pass
        # print(var)
        # cnt += 1
        if len(var.filter.keys()) != 1:
            print('Error: more than one filters for a variant. Exit.', file=sys.stderr)
            print(var)
            exit(1)
        elif var.filter.keys()[0] == 'PASS':
            # Only take 'PASS' variants
            match_allele(var, f_panel, f_fasta)
        input(var)
    # print(cnt)


if __name__ == '__main__':
    args = parse_args()
    # fn_happy_vcf = args.vcf
    # fn_panel_vcf = args.panel
    # fn_summary = args.out

    annotate_vcf(
        fn_vcf=args.vcf,
        fn_panel_vcf=args.panel,
        fn_fasta=args.fasta
        fn_summary=args.out
    )
    exit(0)
    
    print('Warning: the indel precision and F1_score measurements are not identical to official hap.py',
           file=sys.stderr)
    if not args.force_run:
        if os.path.exists(args.out):
            print(f'File "{args.out}" exists. Exit.', file=sys.stderr)
            df = pd.read_csv(args.out, sep='\t')
            print(df)
            exit(0)

    summarize_accuracy(
        fn_happy_vcf=args.vcf,
        fn_summary=args.out
    )


