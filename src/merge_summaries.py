import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-1', '--df1',
        help='Path to df1. [None]'
    )
    parser.add_argument(
        '-2', '--df2',
        help='Path to df2. [None]'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output TSV. ["-"]'
    )
    parser.add_argument(
        '-k1', '--key1', default='CHR,POS,TYPE,REF,ALT,AF',
        help='df1 keys for merging. Fields are separated by commas. ["CHR,POS,TYPE,REF,ALT,AF"]'
    )
    parser.add_argument(
        '-k2', '--key2', default='CHR,POS,TYPE,REF,ALT,AF',
        help='df2 keys for merging. Fields are separated by commas. ["CHR,POS,TYPE,REF,ALT,AF"]'
    )
    args = parser.parse_args()
    return args

def merge_summaries(
    fn_df1, fn_df2, fn_out, 
    key1, key2
    ) -> None:
    df1 = pd.read_csv(fn_df1, sep='\t')
    df2 = pd.read_csv(fn_df2, sep='\t')
    key1 = key1.split(',')
    key2 = key2.split(',')
    df = pd.merge(left=df1, right=df2, how='inner', left_on=key1, right_on=key2)
    df.to_csv(fn_out, sep='\t', index=None) # float_format='%.6f', 


if __name__ == '__main__':
    args = parse_args()
    merge_summaries(
        fn_df1=args.df1, fn_df2=args.df2, fn_out=args.out, 
        key1=args.key1, key2=args.key2
    )
