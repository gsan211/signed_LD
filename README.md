# signed_LD

Collection of scripts and data files associated with Sandler et al. 2021 MBE (Patterns and Causes of Signed Linkage Disequilibria in Flies and Plants)

The "raw" data (annotated VCF files) used for this study are available at https://zenodo.org/record/4895505
The scripts can be used to parse out the data and calculate LD in a few different ways

## site_lists
These two directories contain the lists of sites annotated/used in the paper (based on SNPeff annotations in the VCF files).
Sample SNP eff code is included in the capsella_scripts directory

The SNPeff annotations were parsed with grep to obtain the sites of interest. Have to be careful of multiple annotations when doing this (eg. if you want variants in introns (but not LOF variants), you have to first exlude sites with a "HIGH" predicted effect (this removes splice disrupting variants etc.), then include sites with "intronic_variant" tags.

These site list files are tab delimited and can be used with BCFtools to parse out the sites of interest from the VCFs.

## capsella_scripts
Scripts for calculating LD in Capsella dataset

outlier_samples.txt 

File contains the individuals removed from the VCF for LD analyses 

prepare_vcf_for_R_capsella.sh

Code of how the uploaded VCF was annotated
Code to remove outlier individuals from VCF
Code to filter/parse the VCF (selecting a class of sites from one of the lists in site_lists_capsella directory)


calculate_meanLD_capsella.R

take output file from last step and calculate mean LD for all sites in file, bot site by site and in 100Kb blocks

## drosophila_scripts
Scripts for calculating LD in drosophila dataset

drosophila_inversions_masked.txt

regions that were masked beause of segregatign inversions for some of the supplemental analyses


prepare_vcf_for_R_drosophila.sh

Code to filter/parse the VCF (selecting a class of sites from one of the lists in site_lists_drosophila directory)
Can also be modified to remove inversion regions (use something like BCFtools filter -T ^ ../tab_delimited_regions.txt)


drosophila_edit_genotypes.sh

Takes output from previous step and edits the genotypes to remove missing bases (otherwise LD calcualtions get screwed up)


calculate_meanLD_drosophila.R

Calculates mean LD from file generated by previous script. See annotations in the same capsella script for more details. 


## simulations
admixture_script.slim

Script that was used to run SLiM simulations with admixture (timing of admixture parameter was altered)


calculate_meanLD_simulations.R

Script that takes VCF output of SLiM script and calculates meanLD for deleterious/neutral mutations


## networks
network IDs and genes included in them for both species. See supplementary table 2 from paper for more detail.
Network LD was calculated just like in the calculate_meanLD_capsella.R/calculate_meanLD_drosophila.R
 scripts (in 100Kb windows)

## other scripts included

plink_script.sh

script for generating distribution of LD values using PLINK. Takes as input VCF of sites of interest (can use prepare_vcf_for_R_drosophila.sh script to generate these VCFs)


LD_permutation_test.R / LD_permutation_test_100Kb_windows.R

Permutation test for testing mean LD (SNP by SNP or in 100Kb windows) against 0.
Takes as input output of drosophila_edit_genotypes.sh script (not vcf but the actual table), shuffles genotypes around and generates null distribution of LD. Formatted for Capsella dataset.






