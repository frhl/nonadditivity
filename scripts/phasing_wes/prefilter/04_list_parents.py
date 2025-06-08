#!/usr/bin/env python3

import hail as hl
import argparse
import random

from ukb_utils import hail_init
from ukb_utils import samples

def main(args):

    # parser
    extract_samples = args.extract_samples
    relation = args.relation.upper()
    families = args.families
    seed = args.seed
    use_parents_from_file = args.use_parents_from_file
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/split_parents.log', 'GRCh38')
    ht = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
    overlapping_samples = ht.f0.collect()

    fam = samples.get_fam()

    # ensure that families are in data
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.IID))
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.MAT))
    fam = fam.filter(hl.literal(overlapping_samples).contains(fam.PAT))
   
    # subset to duos or tris
    fam = fam.filter(hl.literal(relation).contains(fam.REL))

    # random family members based on seed:
    if seed and families:
        # Generate a list of random indices without replacement
        random.seed(int(seed))
        num_rows = fam.count() 
        random_indices = random.sample(range(num_rows), min(num_rows,int(families)))
        random_indices_set = hl.literal(set(random_indices))

        # subset pedigree file
        fam = fam.key_by().add_index(name='row_idx')
        fam = fam.transmute(row_idx = hl.int32(fam.row_idx))
        fam = fam.filter(random_indices_set.contains(fam['row_idx']))

        # export new subsetted fam file
        fam = fam.drop(fam.row_idx)
   
    if use_parents_from_file: 
        parents = set(hl.import_table(use_parents_from_file, no_header=True).f0.collect())
        fam = fam.filter(hl.literal(parents).contains(fam.MAT) & hl.literal(parents).contains(fam.PAT)) 

    # export resulting PED file
    fam.export(out_prefix + ".ped")

    # get parents only
    pids = []
    pids.extend(fam.PAT.collect())
    pids.extend(fam.MAT.collect())
    pids = [x for x in pids if x != "0"]

    # filter original file and write out
    ht_p = ht.filter(hl.literal(pids).contains(ht.f0))
    ht_p.f0.export(out_prefix + "_parents.txt")

    # get parents and children
    pids = []
    pids.extend(fam.PAT.collect())
    pids.extend(fam.MAT.collect())
    pids.extend(fam.IID.collect())
    pids = [x for x in pids if x != "0"]

    # export original file with trios
    ht = ht.filter(hl.literal(pids).contains(ht.f0))
    ht.f0.export(out_prefix + "_all.txt")



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--use_parents_from_file', default=None, help='Subset to parents that are already in specified file and get relationships.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--relation', default=None, help='either "trio" or "duo"')
    parser.add_argument('--families', default=None, help='either "trio" or "duo"')
    parser.add_argument('--seed', default=None, help='either "trio" or "duo"')
    args = parser.parse_args()

    main(args)







