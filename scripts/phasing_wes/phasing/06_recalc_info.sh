# author: Frederik Lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

readonly in_dir="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated"
readonly out_dir="/Phasing/PhasingWES/step2_phase_rare_qc/ligated/info"
readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly threads=8

dx mkdir -p ${out_dir}
for CHR in {1..22} X; do
    for anc in "all"; do
      in_file="${in_dir}/UKB_chr${CHR}.exome_array.shapeit5.ligated.bcf"
      samples="${samples_dir}/UKB.wes.qced.${anc}.samples"
      out_file="UKB_chr${CHR}.exome_array.shapeit5.ligated.${anc}.tags.vcf.gz"
      if [[ $(dx_file_exists "${out_dir}/${out_file}") -eq 0 ]]; then
        if [[ "${anc}" == "all" ]]; then
          dx run app-swiss-army-knife -icmd="
             bcftools +fill-tags --threads ${threads} ${in_file} -Ob -o tmp0.bcf -- -t all &&
             bcftools view --drop-genotypes tmp0.bcf -Oz -o ${out_file} &&
             rm tmp*
          " \
          --instance-type mem1_hdd1_v2_x4 \
          --folder=".${out_dir}" \
          --name recalc_info_chr${CHR}_${anc} \
          --priority normal -y
        else
          dx run app-swiss-army-knife -icmd="
             bcftools view ${in_file} --force-samples --samples-file ${samples} -Ob -o tmp0.bcf &&
             bcftools +fill-tags --threads ${threads} tmp0.bcf -Oz -o tmp1.bcf -- -t all &&
             bcftools view --drop-genotypes tmp1.bcf -Oz -o ${out_file} &&
             rm tmp*
          " \
          --instance-type mem1_hdd1_v2_x8 \
          --folder=".${out_dir}" \
          --name recalc_info_${anc}_chr${CHR}_${anc} \
          --depends-on job-GjVKBV8Jg8Jp5P2fbbk2Qvy6 \
          --priority normal -y
        fi
      else
        >&2 echo "${out_file} already exists. Skipping.."
      fi
    done
done

