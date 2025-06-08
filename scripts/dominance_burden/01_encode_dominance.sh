# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=24
readonly out_dir="/wes_ko_ukbb/data/phased/dominance_encoding/"
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/export_alt_qced"
readonly gene_map_dir="/mnt/project/wes_ko_ukbb/data/genesets/dominance_weights"
readonly docker="/docker/call_chets_0.1.14.tar.gz"

dx mkdir -p ${out_dir}
#for pop in "eur" "afr" "eas" "sas"; do
for pop in "eur"; do
  for af in "05"; do
      for CHR in {1..20} 22; do
        gene_map="${gene_map_dir}/UKB.chr${CHR}.phased_sites.spliceai=0.50_cadd=28.1_revel=0.773.canonical.gene_map.txt"
        input_vcf="${in_dir}/UKB.wes.chr${CHR}.phased.qc_final.eur.af${af}.popmax.variants.bcf"
        out_prefix="UKB.wes.chr${CHR}.phased.qc_final.eur.af${af}.popmax.variants.group_dominance_scaling.gene_map"
        dx ls $(wo_mnt_project ${input_vcf})
  #      if [[ $(dx_file_exists "${out_dir}/${out_prefix}.vcf.gz") -eq 0 ]]; then
          dx run app-swiss-army-knife \
            -iimage_file=${docker}\
            -icmd="
              transform --input ${input_vcf} --all-info --gene-map ${gene_map} --set-variant-id --scale-dosages | bgzip > ${out_prefix}.vcf.gz &&
              tabix -C ${out_prefix}.vcf.gz &&
              zcat ${out_prefix}.vcf.gz |  grep -v '##' | cut -f3 | tail -n+2 > ${out_prefix}.sites.txt &&
              echo '$(date)!'
            "\
            --instance-type mem1_ssd1_v2_x2 \
            --folder=".${out_dir}" \
            --priority low \
            --name encode_dominance_c${CHR}_${pop} -y
   #       else
   #         >&2 echo "${out_prefix}.tar.gz already exists. Skipping.."
   #       fi
    done
  done
done


      
