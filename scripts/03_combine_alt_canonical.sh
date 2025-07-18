# author: frederik lassen

set -o errexit
set -o nounset

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
#readonly group="brava_s50"
for pop in "eur" "eas" "afr" "sas"; do
#for pop in "eur"; do
  out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
  dx mkdir -p ${out_dir}
  #for anno in "damaging_missense_lc" "damaging_missense_spliceai_020" "damaging_missense_spliceai_050" "damaging_missense_spliceai_080"; do
  #for anno in "pLoF" "nonsynonymous" "other_missense" "damaging_missense" "pLoF_damaging_missense" "synonymous" "non_coding"; do
  #for anno in "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous"; do
  for anno in "pLoF"; do
    for pp in "0.90"; do
      for af in  "05"; do
        #in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/by_chrom_canonical"
        in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/by_chrom_canonical"
        #in="${in_dir}/UKB.wes.chrCHR.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.txt.gz"
        in="${in_dir}/UKB.wes.chrCHR.phased.full_qcx2.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.txt.gz"
        #out="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.txt.gz"
        out="UKB.wes.merged.phased.full_qcx2.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.txt.gz"
        if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
          in_test=$( wo_mnt_project ${in/CHR/1} )
          if [[ $(dx_is_empty "${in_test}") -eq "0" ]]; then
           dx run app-swiss-army-knife \
             -icmd="
                for idx in {1..22} X; do
                  if [[ \$( zcat ${in/CHR/\$idx} | wc -l) -lt 10 ]]; then
                      echo '${in/CHR/\$idx} does not have any lines!' && 
                      exit 1
                  fi
                done &&
                for idx in {1..22} X; do zcat ${in/CHR/\$idx}; done | sort -k2,3 -V | gzip > ${out} && echo '$(date)'
             " \
             --instance-type mem1_ssd1_v2_x8 \
             --folder=".${out_dir}" \
              --priority high \
              --name combine_chets_${anno}_af${af}_pp${pp} -y
            else
              >&2 echo "Expected ${in_test} to have data! Exiting.." && exit 1
           fi
         else
            >&2 echo "${out} already exists. Skipping."
          fi
        done
     done
  done
done








