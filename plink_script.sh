#script for taking VCF file with sites of interest and using PLINK to estimate signed LD for a samples of those sites

#input VCF into plink and convert into correct format for LD calculations
"/plas1/george.sandler/apps/plink" --make-bed  --vcf "../pop88.syn_max5.vcf" --out ../plink/pop88.syn --double-id  --allow-extra-chr

"/plas1/george.sandler/apps/plink" --bfile ../plink/pop88.syn --recode --tab --out ../plink/pop88.syn --allow-extra-chr

#calculate signed LD. From what I based on some tinkering with mock VCF's, plink assigns positive LD as associations between two minor or two major variants etc.
#this lines up with the way we define LD in the manuscript since we only consider variants below a MAC cut-off of no more than 5
"/plas1/george.sandler/apps/plink" --file ../plink/pop88.syn  --r --out ../plink/pop88.syn --allow-extra-chr --with-freqs
