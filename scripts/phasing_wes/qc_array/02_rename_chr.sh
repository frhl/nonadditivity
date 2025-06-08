# create a file to match chromosome tag between hg19 and hg38 version
#for CHR in {1..22}; do echo "${CHR} chr${CHR}"; done > chr_rename.txt

# upload this file on the UKB RAP
#dx upload chr_rename.txt --path="Phasing/PhasingSNParray/step2_chrrename/"

#rm chr_rename.txt

# rename chromosomes
ANN=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/chr_rename.txt
for CHR in X; do
	VCFF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/full_c$CHR\_b0_v2.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTF --rename-chrs $ANN $VCFF && bcftools index $OUTF" --instance-type mem2_ssd1_v2_x2 --name updatechr_chr$CHR --priority low -y
done

