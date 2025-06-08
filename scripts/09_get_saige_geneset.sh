# author: frederik lassen

# note: this script uses a github repo to get the genests to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="09_get_saige_geneset.R"
readonly rscript_remote="${remote_dir}/09_get_saige_geneset.R"
dx_update_remote ${rscript_remote} ${rscript_local}

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

for pop in "eur"; do
  in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name"
  gene_af_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/gene_af"
  out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}/vcf_plus_plink/force_chr_name/saige_group"
  dx mkdir -p ${out_dir}
  for mode in "dominance" "additive" "recessive"; do
  #for mode in "additive"; do
    #for anno in "nonsynonymous" "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous"; do
    for anno in "pLoF_damaging_missense"; do
       for pp in  "0.90"; do
        for af in "05"; do
          # note: we use recessive (bi-allelic) AFs when mode=Dominance
          gene_ac="${gene_af_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode/dominance/recessive}.maf.txt"
          in_vcf_auto="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.vcf.gz"
          out_prefix_auto="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.chrposrefalt"
          if [[ $(dx_file_exists "${out_dir}/${out_prefix_auto}.vcf.gz") -eq "0" ]]; then 
              dx run app-swiss-army-knife \
                -iimage_file="/docker/rsuite.tar.gz"\
                -icmd="
                   Rscript /mnt/project/${rscript_remote} \
                     --input_ac ${gene_ac} \
                     --canonical_transcripts_only \
                     --out_prefix ${out_prefix_auto} &&
                     echo '!!!!$(date)'
                  "\
                --instance-type mem1_ssd1_v2_x16 \
                --folder=".${out_dir}" \
                --priority normal \
                --name pseudo_variant_vcf_to_saige_group_${mode}_${anno} -y
          else
             >&2 echo "Error! '${in_vcf_auto}' does not exists."
          fi
        done
      done
    done
  done
done





