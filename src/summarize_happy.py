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
        '-s', '--summary',
        help='path to output TSV summary'
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
                print(f'Process {cycle} records', file=sys.stderr)
            cycle += 1

    write_to_tsv(outcome, fn_summary)


if __name__ == '__main__':
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()
    
    print('Warning: the indel precision and F1_score measurements are not identical to official hap.py',
           file=sys.stderr)
    if not args.force_run:
        if os.path.exists(args.summary):
            print(f'File "{args.summary}" exists. Exit.', file=sys.stderr)
            df = pd.read_csv(args.summary, sep='\t')
            print(df)
            exit(0)

    summarize_accuracy(
        fn_happy_vcf=args.vcf,
        fn_summary=args.summary
    )


