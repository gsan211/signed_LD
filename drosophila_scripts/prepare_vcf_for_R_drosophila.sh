bcftoolsg=../bcftools-1.10.2/bin/bcftools

#pull sites of interest from drosophila VCF and write to .txt file as a table
#then input table intro "drosophila_edit_genotypes.sh" script to fix the genotypes
#then LD can be calculated

$bcftoolsg query -f "%ALT\t%CHROM\t%POS[\t%GT]\n" "../ZI_allchrom_clean_bi_samples.ann.vcf" -T ../ZI_utr_max5.list > ../ZI_all_UTR_GTS.txt
