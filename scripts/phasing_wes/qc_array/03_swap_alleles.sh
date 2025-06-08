
# SWAP ALLELES
for CHR in X; do
	VCFF=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/full_c$CHR\_b0_v2.b37.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.swapped.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step3_swapalleles/" -iimage_file="/docker/swaprefalt_0.0.1.tar.gz" -icmd="swapRefAlt_static --input $VCFF --output $OUTF && bcftools index $OUTF" --instance-type mem2_ssd1_v2_x2 --priority low --depends-on job-GYjxyB0Jg8JpGKkv2FZ6y36P --name swapalleles_chr$CHR -y
done

