
readonly out_dir_bmrc="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/vep/dominance_scalings"
readonly in_dir_rap="/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance"


mkdir -p ${out_dir_bmrc}
#for pop in "eur" "afr" "eas" "sas"; do
for pop in "eur"; do
    for CHR in 1; do
      echo "getting scaling for chr${CHR}.."
      file="UKB.wes.chr${CHR}.phased.full_qc.${pop}.af05.popmax.variants.group_dominance_scaling.sites.txt"
      in_file="${in_dir_rap}/${file}"
      out_file="${out_dir_bmrc}/${file}"
      dx cat ${in_file} | sed 's/.*=//' > ${out_file}
    done
done



      
