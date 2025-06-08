# author: Frederik Lassen

set -eu

# Barney has performed filtering of these VCF based on samples.
readonly exome_mnt_in_dir="/mnt/project/Barney/wes/sample_filtered"
readonly array_mnt_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc/support"
readonly support_mnt_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc/support"
readonly merge_mnt_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc/support"
readonly threads=32

# new QCed stuff
readonly qced_variants="/mnt/project/wes_ko_ukbb/data/variants/qced/UKB.wes_qced_variants.combined.txt"

dx mkdir -p ${out_dir}

for CHR in 3 7; do
  array_snps="${array_mnt_dir}/chr${CHR}.array_snps_kept.txt" #
  keep_samples="${merge_mnt_dir}/overlapping_samples.qced.txt"
  bfile="${exome_mnt_in_dir}/ukb_wes_450k.qced.chr${CHR}"
  out_prefix="UKB.chr${CHR}.tmp.full_qc.exome"
  out_expt="${out_dir}/${out_prefix}.bed"
  
  echo ${bfile}
  echo ${array_snps}
  echo ${keep_samples}

  if [[ "$(dx ls ${out_expt} | wc -l)" -eq "0" ]]; then
    dx run app-swiss-army-knife -icmd="
      cat ${qced_variants} | grep 'chr${CHR}:' > keep_variants.txt &&
      plink2 --bed ${bfile}.bed --bim ${bfile}.bim --fam ${bfile}.fam --keep ${keep_samples} --extract keep_variants.txt --exclude ${array_snps} --make-bed --out ${out_prefix} &&
      echo 'o'
    " \
    --instance-type mem2_ssd1_v2_x8 \
    --folder="./${out_dir}" \
    --tag "chr${CHR}" \
    --name filter_qced_wes_chr${CHR} \
    --priority normal -y
  else
    >&2 echo "${out_expt} already exists! Skipping."
  fi
done

