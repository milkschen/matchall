
import argparse
from typing import Tuple
import pysam
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to target VCF. [required]'
    )
    parser.add_argument(
        '-i', '--info-expression', default=None,
        action='append', nargs='+',
        help='Expression for an INFO field, e.g. "AF > 0.5". [None]'
    )
    parser.add_argument(
        '-f', '--format-expression', default=None,
        action='append', nargs='+',
        help='Expression for an FORMAT field, e.g. "BD == FP", "BD[1]=TP". [None]'
    )
    # parser.add_argument(
    #     '--format-idx', default=1, type=int,
    #     help='Index of samples to retrieve format tags. [1]'
    # )
    parser.add_argument(
        '-o', '--out', default='-',
        help='Path to output VCF. Set to "-" to print to stdout. ["-"]'
    )
    parser.add_argument(
        '--happy', action='store_true',
        help='Set for hap.py VCFs. Will only consider variants in confident regions. [False]'
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Set to print debug messages. [False]'
    )
    args = parser.parse_args()
    return args


def parse_expression(exp: str) -> Tuple[str, list, int]:
    '''
    Supported operators (spaces are allowed):
        =, ==, !=, >, >=, <, <=
    e.g. AF > 0.5, BD=FP

    Returns:
        - operator
        - components
        - index
    '''
    re_exp = re.compile(r'[\!\=\<\>\ ]+')
    components = re_exp.split(exp)
    if len(components) != 2:
        raise(ValueError, f'Too many components in expression {exp}. Should be 2.')
    try:
        components[1] = float(components[1])
    except:
        pass

    re_brackets = re.compile('\[(.*?)\]')
    index = re_brackets.findall(components[0])
    if len(index) == 0:
        index = None
    elif len(index) == 1:
        index = int(index[0])
        components[0] = re_brackets.split(components[0])[0]
    else:
        raise(ValueError, f'Invalid component {components[0]}. There should be <= 1 brackets.')

    operator = re_exp.findall(exp)
    if len(operator) != 1:
        raise(ValueError, f'Too many operators in expression {exp}. Should be 1.')
    operator = operator[0].replace(' ', '')
    if operator not in ['=', '==', '!=', '>', '>=', '<', '<=']:
        raise(ValueError, f'Unsupported operator {operator}')
    
    return operator, components, index


def eval_operator(operator, left, right) -> bool:
    if operator in ['=', '==']:
        return left == right
    elif operator == '!=':
        return left != right
    elif operator == '>':
        return left > right
    elif operator == '>=':
        return left >= right
    elif operator == '<':
        return left < right
    elif operator == '<=':
        return left <= right
    else:
        raise(ValueError, f'Unsupported operator {operator}')


def vcf_filter_core(
    var: pysam.VariantRecord, info_exps: list, format_exps: list
    ) -> bool:
    # Check INFO
    if info_exps:
        for opt, (tag, comp), idx in info_exps:
            # If tag is not found, exclude this var
            if not var.info.get(tag):
                if_pass_filter = False
                return False
            elif idx:
                if not eval_operator(operator=opt, left=var.info[tag][idx], right=comp):
                    if_pass_filter = False
                    return False
            # If `idx` is None, loop over all samples
            else:
                for i, _ in enumerate(var.alts):
                    if not eval_operator(operator=opt, left=var.info[tag][i], right=comp):
                        if_pass_filter = False
                        return False
    
    # Check FORMAT
    if format_exps:
        for opt, (tag, comp), idx in format_exps:
            if idx:
                if not var.samples[i].get(tag):
                    if_pass_filter = False
                    return False
                if not eval_operator(operator=opt, left=var.samples[idx][tag], right=comp):
                    if_pass_filter = False
                    return False
            # If `idx` is None, loop over all samples
            else:
                for i, _ in enumerate(var.samples.items()):
                    if not var.samples[i].get(tag):
                        if_pass_filter = False
                        return False
                    if not eval_operator(operator=opt, left=var.samples[i][tag], right=comp):
                        if_pass_filter = False
                        return False
    return True


def vcf_filter(
    fn_vcf: str,
    fn_out: str, 
    info_expression: list,
    format_expression: list,
    happy_vcf: bool=False,
    debug: bool=False
    ) -> None:
    f_vcf = pysam.VariantFile(fn_vcf)
    f_out = pysam.VariantFile(fn_out, 'w', header=f_vcf.header)
    
    info_exps = None
    if info_expression:
        info_exps = [parse_expression(e[0]) for e in info_expression]
    format_exps = None
    if format_expression:
        format_exps = [parse_expression(e[0]) for e in format_expression]

    for var in f_vcf.fetch():
        # In theory we can handle Regions as well, but we haven't handled this
        if happy_vcf:
            if not var.info.get('Regions'):
                continue
            elif 'CONF' not in var.info['Regions']:
                continue
        
        if vcf_filter_core(var=var, info_exps=info_exps, format_exps=format_exps):
            f_out.write(var)
        # if_pass_filter = True
        # # Check INFO
        # if info_exps:
        #     for opt, (tag, comp), idx in info_exps:
        #         # If tag is not found, exclude this var
        #         if not var.info.get(tag):
        #             if_pass_filter = False
        #             break
        #         elif idx:
        #             if not eval_operator(operator=opt, left=var.info[tag][idx], right=comp):
        #                 if_pass_filter = False
        #                 break
        #         # If `idx` is None, loop over all samples
        #         else:
        #             for i, _ in enumerate(var.alts):
        #                 if not eval_operator(operator=opt, left=var.info[tag][i], right=comp):
        #                     if_pass_filter = False
        #                     break
        #     if not if_pass_filter:
        #         continue
        
        # # Check FORMAT
        # if format_exps:
        #     for opt, (tag, comp), idx in format_exps:
        #         if idx:
        #             if not var.samples[i].get(tag):
        #                 if_pass_filter = False
        #                 break
        #             if not eval_operator(operator=opt, left=var.samples[idx][tag], right=comp):
        #                 if_pass_filter = False
        #                 break
        #         # If `idx` is None, loop over all samples
        #         else:
        #             for i, _ in enumerate(var.samples.items()):
        #                 if not var.samples[i].get(tag):
        #                     if_pass_filter = False
        #                     break
        #                 if not eval_operator(operator=opt, left=var.samples[i][tag], right=comp):
        #                     if_pass_filter = False
        #                     break
        #     if not if_pass_filter:
        #         continue
                
        


if __name__ == '__main__':
    args = parse_args()
    vcf_filter(
        fn_vcf=args.vcf, fn_out=args.out, 
        info_expression=args.info_expression, format_expression=args.format_expression,
        happy_vcf=args.happy, debug=args.debug
    )