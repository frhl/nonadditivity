#!/usr/bin/env python3

import hail as hl
import argparse
import random

from ukb_utils import hail_init
from ukb_utils import samples

def main(args):

    # parser
    extract_samples = args.extract_samples
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/split_parents.log', 'GRCh38')
    ht = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
    overlapping_samples = ht.f0.collect()

    fam = samples.get_fam()

    # ensure that families are in data
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.IID))
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.MAT))
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.PAT))
   
    print(fam.REL.collect())
 
    # subset to duos or tris
    #fam = fam.filter(hl.literal(relation).contains(fam.REL))


    # export resulting PED file
    #fam.export(out_prefix + ".ped")



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)

