bcftoolsg=../apps/bcftools-1.10.2/bin/bcftools

#SNPeff command used to generate annotations, requires capsella genome to be added to SNPeff manually (see SNPeff manual)
java -jar "/plas1/george.sandler/apps/snpEff/snpEff.jar" -c /plas1/george.sandler/apps/snpEff/snpEff.config  -v Capsella_rubella "/plas1/george.sandler/capsella/pop88.vcf.gz" > /plas1/george.sandler/capsella/pop88_chrom8.ann.vcf

#command to remove outlier capsella samples from VCF
$bcftoolsg view -Ov "../pop88.ann.vcf" --samples ^143s_clipped,149,165_clipped.bam,27,64_replacement,7s_clipped > "../pop88_sans_outliers.ann.vcf" --force-samples

#pull sites of interest from VCF and create .txt file with site information and matrix of genotypes e.g. synonymous sites
#this is the input for R script: 'calculate_meanLD_capsella.R'
$bcftoolsg filter "../pop88_sans_outliers.ann.vcf" -Ov -e 'TYPE = "indel"' -T "../pop88.miss_syn_max5.list" >"../pop88_SYN.ann.vcf"
$bcftoolsg query "../pop88_SYN.ann.vcf" -f "%CHROM\t%POS[\t%GT]\n" -o /plas1/george.sandler/capsella/fly_retry/pop88_SYN_gts.txt
