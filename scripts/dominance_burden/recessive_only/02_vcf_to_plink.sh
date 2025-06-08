# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=24
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/recessive_encoding/"
readonly out_dir="/wes_ko_ukbb/data/phased/recessive_encoding/"

dx mkdir -p ${out_dir}
for pop in "eur"; do
  for af in "05"; do
      for CHR in {1..22}; do
        input_vcf="${in_dir}/UKB.wes.chr${CHR}.phased.qc_final.eur.af${af}.popmax.variants.recessive.vcf.gz"
        out_prefix="UKB.wes.chr${CHR}.phased.qc_final.eur.af${af}.popmax.variants.recessive"
        dx ls $(wo_mnt_project ${input_vcf})
        if [[ $(dx_file_exists "${out_dir}/${out_prefix}.bed") -eq 0 ]]; then
          dx run app-swiss-army-knife \
            -icmd="
              plink2 --vcf ${input_vcf} dosage=DS --make-bed --out ${out_prefix} &&
              echo '--$(date)' 
            "\
            --instance-type mem1_ssd1_v2_x2 \
            --folder=".${out_dir}" \
            --priority low \
            --name vcf_to_plink_recessive_c${CHR}_${pop} -y
          else
            >&2 echo "${out_prefix}.tar.gz already exists. Skipping.."
          fi
    done
  done
done


      
