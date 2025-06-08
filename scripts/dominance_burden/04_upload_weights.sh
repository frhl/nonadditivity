#!/usr/bin/env bash

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

readonly revel_cutoff="0.773"
readonly cadd_cutoff="28.1"
readonly spliceai_cutoff="0.50"
readonly annotation="spliceai=${spliceai_cutoff}_cadd=${cadd_cutoff}_revel=${revel_cutoff}"

readonly in_dir_bmrc="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/vep/saige_rescaled_weights/af05"
readonly out_dir_rap="/wes_ko_ukbb/data/genesets/dominance_weights/af05"
dx mkdir -p ${out_dir_rap}

for ac_pop in "eur"; do
  #file="UKB.combined.phased_sites.spliceai=0.50_cadd=28.1_revel=0.773.saige_with_${ac_pop}_weights.txt.gz"
  #dx_update_remote "${out_dir_rap}/${file}" "${in_dir_bmrc}/${file}"

  for chr in {1..22} X; do
    file="UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.af05.popmax.variants.${annotation}.min_mac_5.saige_with_${ac_pop}_weights.txt.gz"
    dx_update_remote "${out_dir_rap}/${file}" "${in_dir_bmrc}/${file}"
    file2="UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.af05.popmax.variants.${annotation}.min_mac_5.saige.txt.gz"
    dx_update_remote "${out_dir_rap}/${file}" "${in_dir_bmrc}/${file}"
    
    
    #file="UKB.chr${chr}.phased_sites.${annotation}.saige_with_${ac_pop}_weights.txt.gz"
    #file2="UKB.chr${chr}.phased_sites.${annotation}.saige.txt.gz"
    #dx_update_remote "${out_dir_rap}/${file2}" "${in_dir_bmrc}/${file2}"
    #file3="UKB.chr${chr}.phased_sites.${annotation}.canonical.gene_map.txt"
    #dx_update_remote "${out_dir_rap}/${file3}" "${in_dir_bmrc}/${file3}"
  done
done


