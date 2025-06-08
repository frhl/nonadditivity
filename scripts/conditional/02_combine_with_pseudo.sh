# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly threads=24

# Parameters for pseudo-variants
readonly af="05"
readonly pp="0.90"
readonly pop="eur"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly mode="additive"
readonly in_pseudo_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"

# imputed variants and out
readonly in_imputed_dir="/mnt/project/wes_ko_ukbb/data/conditional/imputed/v2"
readonly out_dir="/wes_ko_ukbb/data/conditional/combined_final/v2"

# parameters for imputed data
readonly maf="01"

dx mkdir -p ${out_dir}
for CHR in {1..22}; do
    for anno in "pLoF" "pLoF_damaging_missense" "synonymous"; do
      pseudo="${in_pseudo_dir}/UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto"
      imputed="${in_imputed_dir}/UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2"
      out_prefix="UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.gel_imputed.maf${maf}.chr${CHR}.v2"
      if [[ $(dx_file_exists "${out_dir}/${out_prefix}.bed") -eq "0" ]]; then
        dx run app-swiss-army-knife \
            -icmd="
               zcat ${pseudo}.vcf.gz | sed -e 's/chr1/chr${CHR}/g' | bgzip > tmp.vcf.gz &&
               tabix -C tmp.vcf.gz &&
               bcftools concat tmp.vcf.gz ${imputed}.bcf -Oz -o ${out_prefix}.vcf.gz &&
               tabix -C ${out_prefix}.vcf.gz &&
               plink2 --vcf ${out_prefix}.vcf.gz dosage=DS --make-bed --out ${out_prefix} &&
               rm ${out_prefix}.vcf.* &&
               rm tmp.vcf*
          "\
            --instance-type mem1_ssd1_v2_x2 \
            --folder=".${out_dir}" \
            --priority low \
            --name c${CHR}_mrg_imputed -y
     else
       >&2 echo "${out_prefix}.bed already exists! Skipping.."
     fi
   done
done

