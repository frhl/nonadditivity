# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly threads=24

# Parameters for pseudo-variants
readonly af="05"
readonly pp="0.90"
readonly pop="eur"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly mode="dominance"
readonly in_pseudo_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"

# imputed variants and out
readonly in_imputed_dir="/mnt/project/wes_ko_ukbb/data/conditional/imputed/recoded/v2"
readonly out_dir="/wes_ko_ukbb/data/conditional/combined_final/v2"

# parameters for imputed data
readonly maf="01" # correspondsto 0.01

dx mkdir -p ${out_dir}
for CHR in {1..22}; do
    for anno in "pLoF" "pLoF_damaging_missense"; do
      recoded_pseudo="${in_pseudo_dir}/UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto"
      recoded_imputed="${in_imputed_dir}/UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2.dominance"
      out_prefix="UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.gel_imputed.maf${maf}.chr${CHR}.v2.dominance"
      dx ls $(wo_mnt_project ${recoded_imputed}.vcf.gz)
      dx ls $(wo_mnt_project ${recoded_pseudo}.vcf.gz)
      if [[ $(dx_file_exists "${out_dir}/${out_prefix}.vcf.gz") -eq "0" ]]; then
        dx run app-swiss-army-knife \
            -icmd="
               zcat ${recoded_pseudo}.vcf.gz | sed -e 's/chr1/chr${CHR}/g' | bgzip > tmp.vcf.gz &&
               tabix -C tmp.vcf.gz &&
               bcftools concat tmp.vcf.gz ${recoded_imputed}.vcf.gz -Oz -o ${out_prefix}.vcf.gz &&
               tabix -C ${out_prefix}.vcf.gz &&
               rm tmp.vcf*
          "\
            --instance-type mem1_ssd1_v2_x2 \
            --folder=".${out_dir}" \
            --priority high \
            --name c${CHR}_combine_with_dom_pseudo_${anno} -y
     else
       >&2 echo "${out_prefix}.bed already exists! Skipping.."
     fi
   done
done

