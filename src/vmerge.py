import matchall
import argparse
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [required]'
    )
    parser.add_argument(
        '-r', '--ref', default=None,
        help='Path to the reference FASTA (FAI index is required). [None]'
    )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output VCF. Set to "-" to print to stdout. ["-"]'
    )
    args = parser.parse_args()
    return args


def vmerge(fn_vcf, fn_out, fn_fasta):
    f_vcf = pysam.VariantFile(fn_vcf)
    f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
    buffer_var = None

    # #TODO
    # if gt not phased:
    #     exit()

    for var in f_vcf.fetch():
        if not buffer_var:
            buffer_var = var
            continue
        # If var overlaps with buffer_var
        if var.start > buffer_var.start + max([len(a) for a in buffer_var.alleles]):
            f_out.write(buffer_var)
            buffer_var = None
        else:
            buffer_var = matchall.vmerge_variant_pair(buffer_var, var)



if __name__ == '__main__':
    # We use python3.6 because we require dict is ordered.
    MIN_PYTHON = (3, 6)
    if sys.version_info < MIN_PYTHON:
        sys.exit('Python %s.%s or later is required.\n' % MIN_PYTHON)

    args = parse_args()
    vmerge(fn_vcf=args.vcf, fn_out=args.out, fn_fasta=args.ref)
