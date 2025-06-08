#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from ukb_utils import hail_init


def main(args):

    # parser
    input_path = args.input_path
    out_prefix = args.out_prefix
    cutoff = float(args.cutoff)

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    ht = hl.read_table(input_path)
    # Restrict to variants in gnomAD with allele frequency > 0.01.
    ht = ht.filter((ht.popmax[0].AF > cutoff) & (hl.len(ht.filters) == 0))
    ht = ht.annotate(
        variant = hl.str(":").join([ht.locus.contig, hl.str(ht.locus.position), ht.alleles[0], ht.alleles[1]])
    )
    ht = ht.key_by()
    ht = ht.select(ht.variant)
    ht.export(out_prefix + ".tsv.bgz", header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--cutoff', default=None, help='cutoff')
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--out_prefix', default=None,help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

