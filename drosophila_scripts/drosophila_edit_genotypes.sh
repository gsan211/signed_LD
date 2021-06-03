#takes in as input txt file (outputted using bcftools from the dorosphila vcf) and converts genotype ID's such that 1=non-reference allele and 0=reference allele
#neccessary as N's are included as 0/1/2 in the original vcf's which messes up downstream calculations of LD or nucleotide diversity

#format must include first column alt alleles, second column chromosome, third column position of SNP, other columns numeric genotypes of all individuals e.g.
#N,A     2L      8203    1       0       0       0       0       0       2       0 
#will convert N alleles to blank spaces, treating them as NA's in further analyses
#converts ref/alt alleles to 0/1 if necessary
#need to provide /tmp directory to write intermediate files
#might also need to sort final file by coordinates for certain downstream analyses

#provide files for input here
for ind in  ../ZI_snpeff_SYN_gts
do

grep -e "N," /plas1/george.sandler/dpgp3/$ind.txt | cut -f4-200 | awk "{gsub("1","n"); print }" OFS="\t"| awk "{gsub("2","1"); print }"  OFS="\t"> ../tmp/N_c_gts_edited.txt
grep -e ",N" /plas1/george.sandler/dpgp3/$ind.txt | cut -f4-200 | awk "{gsub("2","n"); print }"  > ../tmp/c_N_gts_edited.txt

grep -e "N," /plas1/george.sandler/dpgp3/$ind.txt | cut -f1-3 > ../tmp/N_c_names_edited.txt
grep -e ",N" /plas1/george.sandler/dpgp3/$ind.txt | cut -f1-3 > ../tmp/c_N_names_edited.txt

paste /plas1/george.sandler/dpgp3/tmp/N_c_names_edited.txt ..tmp/N_c_gts_edited.txt > ../tmp/N_c_full.txt
paste /plas1/george.sandler/dpgp3/tmp/c_N_names_edited.txt ../tmp/c_N_gts_edited.txt > ../tmp/c_N_full.txt

awk -F "\t" 'length($1)== 1 {print}' ../$ind.txt OFS="\t"> ../tmp/no_edit_needed.txt

cat ./tmp/N_c_full.txt ../tmp/c_N_full.txt ../tmp/no_edit_needed.txt>  ../${ind}_edited.txt

done
