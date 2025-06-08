# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

# note plink can not deal with chrX VCFs:
# Error: Line 13 of --vcf file is for a chrX, chrM, or fully-haploid variant, and
# has a DS field without a companion GT field to clarify whether each DS value is

for pop in "eur" "eas" "afr" "sas"; do
#for pop in "eur"; do
  samples="${samples_dir}/UKB.wes.qced.${pop}${females}.samples" 
  in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"
  out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"
  for mode in "dominance" "additive" "recessive"; do
  #for mode in "dominance"; do
    #for anno in "other_missense"; do
    #for anno in "damaging_missense_lc" "damaging_missense_spliceai_020" "damaging_missense_spliceai_050" "damaging_missense_spliceai_080"; do
    for anno in "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous" "other_missense"; do
      for pp in "0.90"; do
        for af in "05"; do

          out_prefix_auto="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${mode}.auto"
          in_vcf_auto="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${mode}.auto.vcf.gz"
          
          # used for encoding psuedo-variant that have been dominance encoded and scaled across all sites
          #out_prefix_auto="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.scaled"
          #in_vcf_auto="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.scaled.vcf.gz"
          if [[ $(dx_file_exists "${out_dir}/${out_prefix_auto}.bed") -eq "0" ]]; then 
            test_vcf="$(wo_mnt_project ${in_vcf_auto})"
            if [[ $(dx_file_exists ${test_vcf}) -eq "1" ]]; then
              dx run app-swiss-army-knife \
                -icmd="
                  plink2 --vcf ${in_vcf_auto} dosage=DS --keep ${samples} --make-bed --out ${out_prefix_auto} &&  
                  echo '--$(date)'
                  "\
                --instance-type mem1_ssd1_v2_x8 \
                --folder=".${out_dir}" \
                --priority normal \
                --name vcf_to_plink_${mode}_${anno} -y
            else
               >&2 echo "Error! '${in_vcf_auto}' does not exists."
            fi
          fi
        done
      done
    done
  done
done

#plink2 --vcf ${in_vcf_chrx} dosage=DS --keep ${females} --make-bed --out ${out_prefix_chrx} && 
