#dx build -f saige-universal-step-2

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset
set +x
in_dir="/Phasing/PhasingWES/step4_polish_phase/pp_extract"

for chr in {1..22}; do
  
  # setup out/in
  in_dir="/Phasing/PhasingWES/step4_polish_phase/phase_calling/chunks/chr${chr}"
  out_dir="/Phasing/PhasingWES/step4_polish_phase/phase_calling"
  out="UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin_sub_merged.bin"
  dx mkdir -p ${out_dir}

  # make bash expansion for required files 
  files="$(dx ls ${in_dir} | sed 's/^/-i splitted_binary_files=@\//' |  tr "\n" " ")"
  files="${files//"@"/${in_dir}}"
  files="$( echo ${files} | sed 's/,$//g')"
  
  # only run if the previous step has completed
  files_found="$(dx ls ${in_dir} | wc -l)"
  files_expected="42"
  if [[ $(dx_file_exists "${out_dir}/${out}") -eq 0 ]]; then
      if [[ "${files_found}" -eq "${files_expected}" ]]; then
        dx run bin-merger-applet ${files}  --instance-type "mem2_ssd1_v2_x2" --priority high --destination ".${out_dir}" -y --name "c${chr}_bin_merger"
      else
        >&2 echo "Found ${files_found} files but expected ${files_expected}! Cancelling.."
      fi
  else
    >&2 echo "${out} already exists. Skipping.."
  fi
done


