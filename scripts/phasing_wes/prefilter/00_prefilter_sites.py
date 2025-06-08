#!/usr/bin/env python3

import hail as hl
import argparse
import random

from ukb_utils import hail_init
from ukb_utils import samples


def in_target_expr(mt, target_path, padding: int = 0):
    r'''Returns BooleanExpression indicating whether a variant is in the sequencing target

    If locus overlaps with exome target or the padded target, then `in_target`=True
    or `in_target_{padding}bp`=True respectively.

    mt (MatrixTable): MatrixTable to be annotated/filtered
    target_path: Path to xgen target region file provided by UKBB
    padding (int): Base-pair padding to add to target regions, int
    remove (bool): If `remove`=True, this returns a MatrixTable that is filtered
        to rows where `in_target`=True
    '''
    target = hl.import_bed(target_path, reference_genome='GRCh38')
    if padding > 0:
        start = target.interval.start
        end = target.interval.end
        padded_interval = hl.interval(hl.locus(contig=start.contig, pos=start.position-padding),
                                      hl.locus(contig=end.contig, pos=end.position+padding))
        padded_field = f'padded_interval_{padding}bp'
        target = target.annotate(**{padded_field: padded_interval})
        target = target.key_by(padded_field)

    return hl.is_defined(target[mt.locus])

def in_lcr_expr(mt, lcr_path):
    r'''Returns BooleanExpression indicating whether a variant is in a low-complexity region

    If locus is in low-complexity regions, then the expression is True,
    otherwise the expression is False.

    Low-complexity regions are defined using supplementary data provided in
    Li 2014 (https://academic.oup.com/bioinformatics/article/30/20/2843/2422145)

    mt (MatrixTable): MatrixTable to be annotated/filtered
    remove (bool): If `remove`=True, this returns a MatrixTable that is filtered
        to rows where `in_lcr`=False, removing variants located in low-complexity
        regions
    '''
    lcr = hl.import_bed(lcr_path, reference_genome='GRCh38', force_bgz=True)

    return hl.is_defined(lcr[mt.locus])


def main(args):

    input_path = args.input_path
    target_path = args.target_path
    target_padding = int(args.target_padding)
    lcr_path = args.lcr_path
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')  # from zulip
    if "bim" in input_path:
        ht = hl.import_table(input_path, no_header=True)
        ht = ht.rename({'f0': 'contig', 'f1': 'rsid', 'f2': 'cm_position',
                        'f3': 'position', 'f4': 'allele1', 'f5': 'allele2'})
        ht = ht.key_by(locus=hl.locus(ht.contig, hl.int32(
            ht.position), reference_genome='GRCh38'), alleles=[ht.allele1, ht.allele2])
    elif "vcf" in input_path:
        contig_recoding = {
            str(i): f"chr{i}" for i in list(range(1, 23)) + ['X']}
        mt = hl.import_vcf(input_path, reference_genome='GRCh38',
                           contig_recoding=contig_recoding, force_bgz=True)
        ht = mt.rows().select()
    else:
        ht = hl.read_table(input_path)

    # low complexity region padding 
    is_in_lcr = in_lcr_expr(ht, lcr_path)    

    ht = ht.filter(
        ~is_in_lcr
    )

    # filters
    is_in_target = in_target_expr(ht, target_path, target_padding)
    
    ht = ht.filter(
        is_in_target
    ) 
    
    print(ht.count())
    ht.export(out_prefix + ".bim") 
    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='.bim or VCF for input')
    parser.add_argument('--target_path', default=None, help='Path to target')
    parser.add_argument('--target_padding', default=0, help='Padding around target regions')
    parser.add_argument('--lcr_path', default=None, help='Low-complexity regions path')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)







