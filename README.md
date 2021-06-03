# signed_LD

Collection of scripts and data files associated with Sandler et al. 2021 MBE (Patterns and Causes of Signed Linkage Disequilibria in Flies and Plants)

The "raw" data (annotated VCF files) used for this study are available at https://zenodo.org/record/4895505
The scripts can be used to parse out the data and calculate LD in a few different ways

## site_lists
These two directories contain the lists of sites annotated/used in the paper (based on SNPeff annotations in the VCF files).
Sampel SNP eff code is included in the capsella_scripts directory

The SNPeff annotations were parsed with grep to obtain the sites of interest. Have to be careful of multiple annotations when doing this (eg. if you want variants in introns (but not LOF variants), you have to first exlude sites with a "HIGH" prediected effect (this removes splice disrupting variants etc.), then include sites with "intronic_variant" tags.

These site list files are tab delimited and can be given to BCFtools to parse out the sites of interest from the VCFs.

##capsella_scripts
Scripts for calculating LD in Capsella dataset

the outlier_samples.txt 
File contains the individuals removed from the VCF for LD analyses 

prepare_vcf_for_R_capsella.sh
Code of how the uploaded VCF was annotated
Code to remove outlier individuals from VCF
Code to filter/parse the VCF (selecting a class of sites from one of the lists in site_lists_capsella directory)

calculate_meanLD_capsella.R
take output file from last step and calculate mean LD for all sites in file, bot site by site and in 100Kb blocks




