# author: Frederik Lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

in_dir="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated/info"
out_dir="/Phasing/PhasingWES/step2_phase_rare_qc/ligated/info"
threads=8

for anc in "all"; do
  for CHR in {1..22} X; do
      in_file="${in_dir}/UKB_chr${CHR}.exome_array.shapeit5.ligated.${anc}.tags.vcf.gz"
      out_file="UKB_chr${CHR}.exome_array.shapeit5.ligated.${anc}.tags.txt"
      if [[ $(dx_file_exists "${out_dir}/${out_file}.gz") -eq "0" ]]; then
        dx run app-swiss-army-knife -icmd="
           echo -e 'ID\tAC\tAN\tNS\tAF\tMAF\tAC\tAC_Het\tAC_Hom\tAC_Hemi\tHWE\tExcHet' > ${out_file} &&
           bcftools query -f '%CHROM:%POS:%REF:%ALT\t%INFO\n' ${in_file} |  sed -e 's/;/\t/g' | sed 's/[A-z]*=//g' >> ${out_file} &&
           gzip ${out_file}
        " \
        --instance-type mem1_hdd1_v2_x4 \
        --folder=".${out_dir}" \
        --name recalc_info_text_c${CHR} \
        --priority normal -y
      else
        >&2 echo "File already exists: ${out_dir}/${out_file}"
      fi
  done
done
