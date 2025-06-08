set -e  # Exit on error
set -u  # Exit on undefined variables

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# get main script required for conditional analysis
readonly remote_dir="wes_ko_ukbb/scripts"
readonly bashscript_local="_conditional_dom.sh"
readonly bashscript_remote="${remote_dir}/_conditional_dom.sh"
dx_update_remote ${bashscript_remote} ${bashscript_local}

# upload utils too
readonly bashutils_local="_conditional_utils.sh"
readonly bashsutils_remote="${remote_dir}/_conditional_utils.sh"
dx_update_remote ${bashsutils_remote} ${bashutils_local}

# get parameters of the files
readonly anc="eur"
readonly af="05"
readonly pp="0.90"
readonly pop="eur"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly mode="dominance"
readonly maf="01"

# these are the paths. Note that ANNO and CHR is substituted in during the main loop below
# this is the [0,1,2]->[0,0,1] files that can be run quickly to figure out what variants are valid to test
readonly plink_recessive_dir="/mnt/project/wes_ko_ukbb/data/conditional/imputed/recoded/v2"
readonly plink_recessive_prefix="${plink_recessive_dir}/UKB.gel_imputed.sorted.maf${maf}.chrCHR.eur.final.v2.recessive"

# this is the combined imputed and pseudo-variant line, these are files we actually want to test
readonly vcf_dominance_dir="/mnt/project/wes_ko_ukbb/data/conditional/combined_final/v2"
readonly vcf_dominance_prefix="${vcf_dominance_dir}/UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.ANNO.${mode}.auto.gel_imputed.maf${maf}.chrCHR.v2.dominance"

readonly out_dir="/wes_ko_ukbb/data/conditional/output/dominance/v2"
# THis is the path to pseudo variant file generated with 'call_chets.cpp' so that we can link
# individual variants with pseudo-variants (genes)
readonly pseudo_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
readonly pseudo_file="${pseudo_dir}/UKB.wes.merged.phased.qc_final.${pop}.af${af}.popmax.pp${pp}.${group}.ANNO.txt.gz"

dx mkdir -p ${out_dir}

# Set threshold for (dominance) p-value filtering
readonly p_dom_threshold=0.001
readonly p_add_threshold=0.5

# what is the threshold for considering a nearby marker in the conditional analysis?
readonly cond_p_threshold=0.000001
readonly max_iterations=10

# padding around target gene. Note that this should depend on 
# the imputed data regions that were created in the earlier steps!
readonly padding=500000

# get the markers
readonly cond_marker_file="/wes_ko_ukbb/data/conditional/df2_add_dom_cond_intervals_05Feb2025.txt"
readonly cond_marker_file_local="$(basename ${cond_marker_file})"

echo "${cond_marker_file_local}"
if [[ ! -f ${cond_marker_file_local} ]]; then
  dx cat ${cond_marker_file} > ${cond_marker_file_local}
  echo "${cond_marker_file_local}... downloaded!"
else
  echo "${cond_marker_file_local}... ok!"
fi

while IFS=$'\t' read -r \
    ensembl_transcript_id \
    dataset \
    annotation \
    ensembl_gene_id \
    hgnc_symbol \
    chromosome \
    trait \
    title \
    p_value_add \
    p_value_dom \
    start_position \
    end_position; do
    
    # skip header
    if [[ "${trait}" == "trait" ]]; then
      continue
    fi
    
    # skip header
    #if [[ "${chromosome}" != "1" ]]; then
    #  continue
    #fi

    #echo "Only hand grip strength!"
    #if [[ ${trait} != "FOLR3" ]]; then
    #  continue
    #fi

    #echo "Only hand grip strength!"
    #if [[ ${title} != "hand_grip_strength" ]]; then
    #  continue
    #fi

    # Convert scientific notation to decimal before comparison
    p_value_dom_decimal=$(echo "$p_value_dom" | awk '{printf "%.10f\n", $1}')
    if (( $(echo "$p_value_dom_decimal > $p_dom_threshold" | bc -l) )); then
        continue
    fi
    
    p_value_add_decimal=$(echo "$p_value_add" | awk '{printf "%.10f\n", $1}')
    if (( $(echo "$p_value_add_decimal > $p_add_threshold" | bc -l) )); then
        continue
    fi
    
    echo "Processing $hgnc_symbol - $trait (p_value = $p_value_add, p_value_dom = $p_value_dom)"
    
    # setup path to corresponding file required for conditional analysis
    # and use substitution to get right chrom and/or annotation
    plink_recessive_chr=${plink_recessive_prefix/CHR/${chromosome}}
    plink_recessive_chr_anno=${plink_recessive_chr/ANNO/${annotation}}
    vcf_dominance_chr=${vcf_dominance_prefix/CHR/${chromosome}}
    vcf_dominance_chr_anno=${vcf_dominance_chr/ANNO/${annotation}}
    pseudo_file_anno=${pseudo_file/ANNO/${annotation}}

		#echo "checking..
    #echo ${plink_recessive_chr_anno}"
    #dx ls $(wo_mnt_project ${vcf_dominance_chr_anno}.vcf.gz )
    #dx ls $(wo_mnt_project ${plink_recessive_chr_anno}.bim)

    if [[ $( echo ${dataset} | grep olink | wc -l ) -eq 1 ]]; then
      pheno="${trait}"
      step0_dir="/mnt/project/wes_ko_ukbb/data/saige/step0/vr_5k_olink"
      step1_dir="/mnt/project/wes_ko_ukbb/data/saige/step1/quantitative/olink"
    else
      pheno="${trait}.${title}"
      step0_dir="/mnt/project/wes_ko_ukbb/data/saige/step0/vr_20k"
      step1_dir="/mnt/project/wes_ko_ukbb/data/saige/step1/2df"
    fi

    # deal with hand_grip_strength_right (edge case)
    pheno="${pheno/hand_grip_strength/hand_grip_strength_right}"
    pheno=${pheno/p4079.diastolic_bp/participant.p4079_i0_avg_a0_a1}
    pheno=${pheno/p4080.systolic_bp/participant.p4080_i0_avg_a0_a1}
    pheno=${pheno/fev1_fvc_ratio.FEV1_FVC/fev1_fvc_ratio}
    
    
    # file name will be too long for this trait otherwise
    pheno_short=${pheno}
    pheno=${pheno/p20150.FEV1/p20150.forced_expiratory_volume_in_1_second_fev1_best_measure}

    # setup step 0-1 files
    model_file="${step1_dir}/${pheno}_${anc}.rda"
    var_ratio_file="${step1_dir}/${pheno}_${anc}.varianceRatio.txt"
    grm="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
    grm_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
    
    # Define plink file and output prefix
    out_prefix="nearby_common_conditional.maf_cutoff=${maf}.pop=${anc}.trait=${pheno_short}.symbol=${hgnc_symbol}.transcript=${ensembl_transcript_id}.chrom=${chromosome}.group=${group}.anno=${annotation}.mode=${mode}"
    out="${out_dir}/${out_prefix}.txt.gz"

    if [[ $(dx_file_exists $(wo_mnt_project "${plink_recessive_chr_anno}.bed")) -eq "1" ]]; then
      if [[ $(dx_file_exists "${out}") -eq "0" ]]; then
        dx run app-swiss-army-knife \
          -icmd="
            bash /mnt/project/${bashscript_remote} \
              --pheno '${pheno}' \
              --model_file '${model_file}' \
              --var_ratio_file '${var_ratio_file}' \
              --grm '${grm}' \
              --grm_samples '${grm_samples}' \
              --plink_recessive_file '${plink_recessive_chr_anno}' \
              --vcf_dominance_file '${vcf_dominance_chr_anno}' \
							--out_prefix '${out_prefix}' \
              --start_pos '${start_position}' \
              --end_pos '${end_position}' \
              --chromosome '${chromosome}' \
              --transcript_id '${ensembl_transcript_id}' \
              --max_iterations '${max_iterations}' \
              --threshold '${cond_p_threshold}' \
              --pseudovariant_to_variant_file '${pseudo_file_anno}' \
              --padding '${padding}' &&
              echo '$(date)'
          " \
          --instance-type mem3_ssd1_v2_x4 \
          --folder=".${out_dir}" \
          --priority high \
          --name "cond_dom_${pheno}_${hgnc_symbol}" -y
          #echo "exiting.."
          #exit 1
      else
        >&2 echo "${out} already exists. Skipping."
      fi
      
      echo "Submitted analysis for ${pheno} at $(date)"
    else 
      >&2 echo "${plink_recessive_chr_anno} does not exists. Skipping.."
    fi

done < ${cond_marker_file_local}

echo "All conditional analyses completed at $(date)"
