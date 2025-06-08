# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# Note: an error here without any output can indicate that input file is empty!

set -o errexit
set -o nounset

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly docker="/docker/call_chets_0.2.1.tar.gz"

for pop in "eur" "afr" "eas" "sas"; do
#for pop in "eur"; do
 in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
 out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"
 dx mkdir -p ${out_dir}
 for mode in "additive" "recessive" "dominance"; do
 #for mode in "dominance" "recessive"; do
   #for anno in "damaging_missense_lc" "damaging_missense_spliceai_020" "damaging_missense_spliceai_050" "damaging_missense_spliceai_080"; do
   for anno in "nonsynonymous" "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous" "nonsynonymous" "other_missense"; do
   #for anno in "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous" "nonsynonymous"; do
      for pp in "0.90"; do
        for af in "05"; do
          # ensure that geno is present
          geno="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.txt.gz"
          samples="${samples_dir}/UKB.wes.qced.${pop}.samples"
          
          # submit combined (except for chrX which is females only) 
          females="${samples_dir}/UKB.wes.qced.${pop}.females.samples"
          out_prefix_auto="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${mode}.auto"
          out_prefix_chrx="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${mode}.chrx"
          out="${out_prefix_auto}.vcf.gz"
         
          #dx rm -f "${out_dir}/${out_prefix_auto}.vcf.gz*"
          #dx rm -f "${out_dir}/${out_prefix_chrx}.vcf.gz*"
          #echo "GENO = ${geno} .. File size = $(dx_size_in_bytes $(wo_mnt_project ${geno}))"
          #echo "SAMPLES = ${samples}"
          #echo "OUT = ${out}"
          
          #if [[ 1 -eq 2 ]]; then
          #if [[ $(dx_is_empty $( wo_mnt_project ${geno}) ) -eq "0" ]]; then
          #  if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
              echo "Submitting regular encoding.."
              dx run app-swiss-army-knife \
                -iimage_file=${docker}\
                -icmd="
                  set -x &&
                  zcat ${geno} | grep -wE 'chr[0-9]+' > chr1_22.txt.gz &&
                  zcat ${geno} | grep -w  'chrX' > chrX.txt.gz &&
                  encode_vcf --input chr1_22.txt.gz --min-ac 1 --force-chr-out-name chr1 --samples ${samples} --mode ${mode} --all-info | bgzip > ${out_prefix_auto}.vcf.gz &&
                  encode_vcf --input chrX.txt.gz --min-ac 1 --force-chr-out-name chr1 --samples ${females} --mode ${mode} --all-info | bgzip > ${out_prefix_chrx}.vcf.gz &&
                  tabix -C ${out_prefix_auto}.vcf.gz &&
                  tabix -C ${out_prefix_chrx}.vcf.gz &&
                  rm -f chr1_22.txt.gz &&
                  rm -f chrX.txt.gz &&
                  echo '$(date)' 
                  "\
                --instance-type mem3_ssd1_v2_x8 \
                --folder=".${out_dir}" \
                --priority normal \
                --name encode_${mode}_${anno} -y
           #  else
           #    >&2 echo "${out} already exists. Skipping."
           #  fi
           #else
           # >&2 echo "File is empty: '${geno}' Exiting loop."  && exit 1
           #fi

          out_scaled="${out_prefix_auto}.scaled.vcf.gz"
         
           if [[ "${mode}" == "dominance" ]]; then
           #  if [[ $(dx_is_empty $( wo_mnt_project ${geno}) ) -eq "0" ]]; then
           #   if [[ $(dx_file_exists "${out_dir}/${out_scaled}") -eq "0" ]]; then
                echo "Submitting fully scaled dominance encoding.."
                dx run app-swiss-army-knife \
                  -iimage_file=${docker}\
                  -icmd="
                    set -x &&
                    zcat ${geno} | grep -wE 'chr[0-9]+' > chr1_22.txt.gz &&
                    zcat ${geno} | grep -w  'chrX' > chrX.txt.gz &&
                    encode_vcf --input chr1_22.txt.gz --min-ac 1 --force-chr-out-name chr1 --samples ${samples} --mode ${mode} --global-dom-dosage --all-info | bgzip > ${out_prefix_auto}.scaled.vcf.gz &&
                    encode_vcf --input chrX.txt.gz --min-ac 1 --force-chr-out-name chr1 --samples ${females} --mode ${mode} --global-dom-dosage --all-info | bgzip > ${out_prefix_chrx}.scaled.vcf.gz &&
                    tabix -C ${out_prefix_auto}.scaled.vcf.gz &&
                    tabix -C ${out_prefix_chrx}.scaled.vcf.gz &&
                    rm -f chr1_22.txt.gz &&
                    rm -f chrX.txt.gz &&
                    echo '$(date)'  
                  "\
                  --instance-type mem3_ssd1_v2_x8 \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name encode_scaled_${mode}_${anno} -y
          #     else
          #       >&2 echo "${out} already exists. Skipping."
          #     fi
          #   else
          #    >&2 echo "File is empty: '${geno}' Exiting loop."  && exit 1
          #   fi
          fi
 
        done
      done
    done
  done
done


