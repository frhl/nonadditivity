# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=24
readonly out_dir="/wes_ko_ukbb/data/conditional/imputed/v2"
readonly in_dir="/mnt/project/Bulk/Imputation/Imputation from genotype (GEL)"

# subst to exome capture regions +/- padding
readonly padding="500000" # in bp
readonly input_xgen="/mnt/project/wes_ko_ukbb/data/conditional/df2_add_dom_cond_intervals_05Feb2025.bed"
#readonly input_xgen="/mnt/project/wes_ko_ukbb/xgen_plus_spikein.b38.chr_prefix.bed"

dx ls $(wo_mnt_project ${input_xgen} )

# NOTE mem2_ssd1_v2_x2 has too low disk space for chr3!

dx mkdir -p ${out_dir}
for CHR in {1..22}; do
  for pop in "eur"; do
    samples_qc="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.${pop}.samples"
    input_bgen="${in_dir}/ukb21008_c${CHR}_b0_v1.bgen"
		samples="${in_dir}/ukb21008_c${CHR}_b0_v1.sample"
    out_prefix="UKB.gel_imputed.sorted.maf01.chr${CHR}.${pop}.final.v2"
    out="${out_prefix}.bcf"
    if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
      dx run app-swiss-army-knife \
        -icmd="
					awk '{print \$1, \$1}' \"${samples_qc}\" > samples.txt &&
					plink2 --bgen \"${input_bgen}\" ref-first  \
           --indiv-sort n \
           --extract bed1 ${input_xgen} \
           --bed-border-bp ${padding} \
           --sample \"${samples}\" \
           --keep samples.txt \
           --geno 0.01 \
           --maf 0.01 \
           --make-bed \
           --out ${out_prefix} &&
				  rm samples.txt &&
          plink2 --bfile ${out_prefix} \
           --output-chr chrM \
           --export bcf vcf-dosage=DS-only id-paste=iid \
           --out ${out_prefix} &&
         rm ${out_prefix}.bed &&
         rm ${out_prefix}.bim &&
         rm ${out_prefix}.fam &&
         rm ${out_prefix}.log
            "\
        --instance-type mem3_ssd1_v2_x4 \
        --folder=".${out_dir}" \
        --priority normal \
        --name c${CHR}_${pop}_imputed_qc -y
    else
      >&2 echo "${out} already exists!"
    fi
  done
done
#   
