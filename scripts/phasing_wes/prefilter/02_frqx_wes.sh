# author: Frederik Lassen

set -eu

readonly qced_eur="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.samples"
readonly exome_mnt_in_dir="/mnt/project/Barney/wes/sample_filtered"
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc/support"
readonly threads=24

dx mkdir -p ${out_dir}

for CHR in 22; do
  bfile="${exome_mnt_in_dir}/ukb_wes_450k.qced.chr${CHR}"
  out_prefix="UKB.chr${CHR}.full_qc.exome"
    dx run app-swiss-army-knife -icmd="
      plink --bfile ${bfile} --freqx --keep  <( cat ${qced_eur} | awk '{print \$0 \"\t\" \$0}') --out ${out_prefix}.eur &&
      gzip ${out_prefix}.eur.frqx &&
      rm *.log &&
      rm *.nosex
    " \
    --instance-type mem2_ssd1_v2_x8 \
    --folder="./${out_dir}" \
    --tag "chr${CHR}" \
    --name frqx_wes_chr${CHR} \
    --priority normal -y
done

