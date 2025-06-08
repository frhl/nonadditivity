# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# Note: an error here without any output can indicate that input file is empty!

set -o errexit
set -o nounset

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly docker="/docker/call_chets_0.1.15.tar.gz"

# we only want to recode additive pseudo-variant!
readonly mode="additive" 

#for pop in "eur" "afr" "eas" "sas"; do
for pop in "eur"; do
 in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/vcf_plus_plink"
 out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/vcf_plus_plink/full_scaled_dominance"
 dx mkdir -p ${out_dir}
   for anno in "pLoF_damaging_missense"; do
   #for anno in "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous" "nonsynonymous"; do
      for pp in "0.90" ; do
        for af in "05"; do
          input_vcf="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.vcf.gz"
          out_prefix="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.auto.fully_scaled_dominance"
          out="${out_prefix}.vcf.gz"
          if [[ $(dx_is_empty $( wo_mnt_project ${input_vcf}) ) -eq "0" ]]; then
            if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
              dx run app-swiss-army-knife \
                -iimage_file=${docker}\
                -icmd="
                    zcat ${input_vcf} | head &&
                    zcat ${input_vcf} | wc -l &&
                    transform --input ${input_vcf} --scale-dosages | bgzip > ${out_prefix}.vcf.gz &&
                    tabix -C ${out_prefix}.vcf.gz &&
                    zcat ${out_prefix}.vcf.gz |  grep -v '##' | cut -f3 | tail -n+2 > ${out_prefix}.sites.txt &&
                    echo '$(date)!' 
                  "\
                --instance-type mem1_ssd1_v2_x8 \
                --folder=".${out_dir}" \
                --priority normal \
                --name scale_pseudo_variants_${anno} -y
             else
               >&2 echo "${out} already exists. Skipping."
             fi
           else
            >&2 echo "File is empty: '${input_vcf}' Exiting loop."  && exit 1
           fi

         done
      done
  done
done
