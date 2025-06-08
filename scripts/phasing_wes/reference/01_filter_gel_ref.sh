# author: frederik lassen

set -eu

threads=16

# get dir to reference and array SNPs we carry forward
ref_dir="/mnt/project/Bulk/Imputation/Imputation' 'from' 'genotype' '\(GEL\)"
array_dir="/mnt/project/Phasing/PhasingWES/step0_merge/support"
out_dir="/Phasing/PhasingWES/reference_panel"

# subset to the same guys going into phasing
keep_samples="/mnt/project/Phasing/PhasingWES/step0_merge/overlapping_samples.qced.txt"
exclude_samples="/mnt/project/Phasing/PhasingWES/step0_merge/ukb_wes_200k_trios_parents.txt"
memory="$((60*1000))"
dx mkdir -p ${out_dir}

for CHR in 18; do
  ref_bgen="${ref_dir}/ukb21008_c${CHR}_b0_v1.bgen" 
  ref_sam="${ref_dir}/ukb21008_c${CHR}_b0_v1.sample"  
  array_snps="${array_dir}/chr${CHR}.array_snps_kept.txt"
  out_prefix="UKB.array.phased.genotypes.from.gel.chr${CHR}"
  dx run app-swiss-army-knife -icmd="
    cat ${array_snps} | cut -d':' -f1,2 > pos.txt &&
    qctool -g ${ref_bgen} -s ${ref_sam} -filetype bgen -incl-samples ${keep_samples} -excl-samples ${exclude_samples} -incl-positions pos.txt -ofiletype vcf -og ${out_prefix}.vcf.gz &&
    rm pos.txt
  "\
    --instance-type mem2_ssd1_v2_x16\
    --folder=".${out_dir}"\
    --name filter_gel_ref_c${CHR}\
    --priority normal -y
done

 #   plink2 --bfile tmp --extract ${array_snps} --keep ${keep_samples} --remove ${exclude_samples} --output-chr chrMT --recode vcf bgz --out ${out_prefix} &&
    #rm tmp* &&
    #plink2 --bgen ${ref_bgen} ref-first --sample ${ref_sam} --memory ${memory} --keep ${keep_samples} --remove ${exclude_samples} --debug  --make-pgen --out ${out_prefix}
   #plink2 --bgen ${ref_bgen} ref-first --sample ${ref_sam} --memory ${memory} --extract ${array_snps} --keep ${keep_samples} --remove ${exclude_samples} --make-pgen --out ${out_prefix}
   #plink2 --bgen ${ref_bgen} ref-first --sample ${ref_sam} --extract ${array_snps} --keep ${keep_samples} --remove ${exclude_samples} --output-chr chrMT --recode vcf bgz --out ${out_prefix}
 #plink --bgen ${ref_bgen} ref-unknown --sample ${ref_sam} --memory ${memory} --extract ${array_snps} --keep ${keep_samples} --remove ${exclude_samples} --output-chr chrMT --recode vcf bgz --out ${out_prefix}




