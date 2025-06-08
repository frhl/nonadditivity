# create local folders
dx mkdir -p Phasing/PhasingSNParray/step1_dataqc/
dx mkdir -p Phasing/PhasingSNParray/step2_chrrename/
dx mkdir -p Phasing/PhasingSNParray/step3_swapalleles/
dx mkdir -p Phasing/PhasingSNParray/step4_liftover/
dx mkdir -p Phasing/PhasingSNParray/step5_phasing/

#Get SNP list
#dx run app-swiss-army-knife -iin="/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="cat ukb_snp_qc.txt | awk '{ print \$1, \$159; }' > SNPlist.unfiltered.txt && cat SNPlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1; }' > SNPlist.filtered.QC.txt" --instance-type mem1_ssd1_v2_x2 --name qc_snp --priority normal -y

#Get sample lists
#dx run app-swiss-army-knife --folder="/Phasing/PhasingSNParray/step1_dataqc/" -iin="/Bulk/Genotype Results/Genotype calls/ukb_sqc_v2.txt" -iin="/Bulk/Genotype Results/Genotype calls/ukb22418_c1_b0_v2.fam" -icmd="cat ukb_sqc_v2.txt | cut -d ' ' -f 66 > in_phasing.txt && cat ukb22418_c1_b0_v2.fam | cut -d ' ' -f 1 > samples.tmp.txt && paste samples.tmp.txt in_phasing.txt > INDlist.unfiltered.txt && cat INDlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1, \$1; }' > INDlist.filtered.QC.txt && rm samples.tmp.txt && rm in_phasing.txt" --instance-type mem1_ssd1_v2_x4 --priority normal --name qc_sample -y

#Filter each chromosome
for CHR in X; do
	dx run app-swiss-army-knife -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bed" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bim" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.fam" -iin="/Phasing/PhasingSNParray/step1_dataqc/INDlist.filtered.QC.txt" -iin="/Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep INDlist.filtered.QC.txt --extract SNPlist.filtered.QC.txt --export vcf bgz --out full_c${CHR}_b0_v2 && bcftools index full_c${CHR}_b0_v2.vcf.gz" --instance-type mem1_ssd1_v2_x2 --name qc_chr${CHR} --priority low -y
done



