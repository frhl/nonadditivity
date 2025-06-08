# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

set -eu

readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/export_alt_qced"
readonly docker="/docker/call_chets_0.1.9.tar.gz"

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly vep_dir="/mnt/project/wes_ko_ukbb/data/vep_loftee/csqs_split_annotation/${group}"
readonly instance_type="mem1_ssd1_v2_x16"

#for pop in "afr" "eas" "sas"; do
for pop in "eur"; do
  out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/by_chrom_canonical"
  dx mkdir -p ${out_dir}
  for anno in "damaging_missense_lc" "damaging_missense_spliceai_020" "damaging_missense_spliceai_050" "damaging_missense_spliceai_080"; do
    for pp in "0.50" "0.90"; do
      for af in "05"; do
        for CHR in {1..22} X; do
          
          # note: this is by all transcripts - we filter later on
          geno="${in_dir}/UKB.wes.chr${CHR}.phased.full_qc.${pop}.af${af}.popmax.txt.gz"
          vep="${vep_dir}/UKB.chr${CHR}.phased_sites.${group}.split.canonical.txt.gz"
          out="UKB.wes.chr${CHR}.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.txt.gz"
          
          echo "GENO = ${geno}"
          echo "VEP =  ${vep}"
          echo "OUT = ${out_dir}/${out}"
          echo " -- File size = $(dx_size_in_bytes ${out_dir}/${out})"

          if [[ "$(dx_file_exists ${out_dir}/${out})" -eq 1  ]]; then
            echo "Output was not exported correctly. Removing and re-submitting!"
            dx rm ${out_dir}/${out}
          fi
          
          #if [[ "$(dx_file_exists ${out_dir}/${out})" -eq 1  && "$(dx_size_in_bytes ${out_dir}/${out} )" -lt 200 ]]; then
          #  echo "Output was not exported correctly. Removing and re-submitting!"
          #  dx rm ${out_dir}/${out}
          #fi
          
          if [[ "$(dx_is_empty $( wo_mnt_project ${geno} ))" -eq 0 ]]; then
            if [[ $(dx_file_exists "${out_dir}/${out}") -eq 0 ]]; then
              if [[ "${anno}" == "nonsynonymous" ]]; then
                dx run app-swiss-army-knife \
                  -iimage_file=${docker}\
                  -icmd="
                      zcat ${vep} | grep -Ew '(pLoF)|(damaging_missense)|(other_missense)' | cut -f1,7 | gzip > map.txt.gz &&
                      zcat ${geno} | awk '\$4>=${pp} || \$4==\"\" || \$4==\".\"' | gzip > geno.txt.gz &&
                      call_chets --geno geno.txt.gz --gene-map map.txt.gz --show-variants | gzip > ${out}
                      rm map.txt.gz &&
                      rm geno.txt.gz &&
                      echo 'xp1zx ${anno}'
                  " \
                  --instance-type ${instance_type} \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name annotate_${anno}_chets_chr${CHR} -y
              elif [[ "${anno}" == "pLoF_damaging_missense" ]]; then
                dx run app-swiss-army-knife \
                  -iimage_file=${docker}\
                  -icmd="
                      zcat ${vep} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f1,7 | gzip > map.txt.gz &&
                      zcat ${geno} | awk '\$4>=${pp} || \$4==\"\" || \$4==\".\"' | gzip > geno.txt.gz &&
                      call_chets --geno geno.txt.gz --gene-map map.txt.gz --show-variants | gzip > ${out}
                      rm map.txt.gz &&
                      rm geno.txt.gz &&
                      echo 'xp1zx ${anno}'
                  " \
                  --instance-type ${instance_type} \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name annotate_${anno}_chets_chr${CHR} -y
              else
               dx run app-swiss-army-knife \
                  -iimage_file=${docker}\
                  -icmd="
                      zcat ${vep} | grep -Ew ${anno} | cut -f1,7 | gzip > map.txt.gz &&
                      zcat ${geno} | awk '\$4>=${pp} || \$4==\"\" || \$4==\".\"' | gzip > geno.txt.gz &&
                      call_chets --geno geno.txt.gz --gene-map map.txt.gz --show-variants | gzip > ${out}
                      rm map.txt.gz &&
                      rm geno.txt.gz &&
                      echo 'okz2af ${anno}'
                  " \
                  --instance-type ${instance_type} \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name c${CHR}_${pop}_encode_${anno} -y
                fi
              else
                >&2 echo "${out} already exists! Skipping.."
              fi
            else
              >&2 echo "${geno} does not exists or is empty!"
            fi

          done
        done
     done
  done
done




