# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

#for pop in "eur" "afr" "eas" "sas"; do
for pop in "eur"; do
 in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/vcf_plus_plink"
 out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced/${pop}/${group}/gene_af"
 dx mkdir -p ${out_dir}
 for mode in "recessive" "additive"; do
   #for anno in "pLoF"; do
   for anno in "pLoF" "pLoF_damaging_missense" "damaging_missense" "synonymous" "nonsynonymous" "other_missense"; do
      for pp in "0.50" "0.90" ; do
        for af in "05"; do
          input_vcf_auto="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto.vcf.gz"
          input_vcf_chrx="${in_dir}/UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.chrx.vcf.gz"
          out_prefix="UKB.wes.merged.phased.full_qc.${pop}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.maf"
          out="${out_prefix}.txt" 
         
#          if [[ $(dx_is_empty $( wo_mnt_project ${input_vcf_auto}) ) -eq "0" ]]; then
#            if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
              dx run app-swiss-army-knife \
                -icmd="
                    bcftools query -f \"%CHROM:%POS:%REF:%ALT;%ID;%INFO\n\" ${input_vcf_auto} >> vcf_content.txt &&
                    bcftools query -f \"%CHROM:%POS:%REF:%ALT;%ID;%INFO\n\" ${input_vcf_chrx} >> vcf_content.txt &&
                    cat vcf_content.txt | cut -d';' -f1-4 | sed -e 's/AC\=//g' | sed -e 's/AN\=//g' > processed.txt &&
                    cat processed.txt | while IFS=\";\" read -r ID GENE AC AN; do
                        AF=\$(awk \"BEGIN {print \$AC / \$AN}\") &&
                        MAF=\$(awk \"BEGIN {if (\$AF < 1 - \$AF) print \$AF; else print 1 - \$AF}\") &&
                        echo -e \"\$ID\t\$GENE\t\$AC\t\$AN\t\$AF\t\$MAF\" >> ${out_prefix}.txt
                    done &&
                    sed -i '1i ID\tGENE\tAC\tAN\tAF\tMAF' processed.txt &&
                    rm processed.txt &&
                    rm vcf_content.txt
                    "\
                --instance-type mem2_ssd1_v2_x4 \
                --folder=".${out_dir}" \
                --priority normal \
                --name gene_af_${anno} -y
  #           else
  #             >&2 echo "${out} already exists. Skipping."
  #           fi
  #         else
  #          >&2 echo "File is empty: '${geno}' Exiting loop."  && exit 1
  #         fi

         done
      done
    done
  done
done


