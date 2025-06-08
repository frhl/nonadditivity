# author: frederik lassen

readonly CHR=21
readonly exome_qualifying_chunk="/mnt/project/Bulk/Exome\ sequences_Alternative\ exome\ processing/Exome\ variant\ call\ files\ \(gnomAD\)\ \(VCFs\)/ukb24068_c${CHR}_b1_v1.vcf.gz"	
readonly array="/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz"
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc"
readonly overlapping_samples="chr${CHR}.overlapping_samples"
readonly qced_samples="/mnt/project/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr21.fam"
readonly extra_qced_samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.combined"

dx mkdir -p ${out_dir}
dx run app-swiss-army-knife -icmd="
  bcftools query -l ${exome_qualifying_chunk} > exome.samples.txt &&
  bcftools query -l ${array} | cut -d '_' -f 1 > array.samples.txt &&
  cat exome.samples.txt array.samples.txt > samples.txt &&
  sort samples.txt | uniq -c | gawk '\$1==2{print \$2}' > overlapping_samples.txt  &&
  rm samples.txt &&
  cat ${qced_samples} | cut -f2 > qced.samples.txt &&
  cat overlapping_samples.txt qced.samples.txt > samples.txt &&
  sort samples.txt | uniq -c | gawk '\$1==2{print \$2}' > overlapping_samples.qced.txt  &&
  rm samples.txt &&
  comm -12 <(sort ${overlapping_samples}) <(sort ${extra_qced_samples}) > overlapping_samples.qced.qced.txt
" \
  --instance-type mem1_ssd1_v2_x4\
  --folder=".${out_dir}"\
  --name WESarray_overlapping_samples\
  --priority normal -y



