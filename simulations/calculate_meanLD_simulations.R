#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

library(matrixStats)
library(data.table)
library(magrittr)
library(windowscanr)
library(dplyr)
library(gdsfmt)
library(SNPRelate)

#read in a file with a list of paths for the vcf's outputed by slim
ll = fread(  "../list_ad1300k.list", header = F)

analyze  <- function(line_filename){
  line = nth(ll$V1, line_filename)
  
#perform PCA on sampled genotypes with SNPRelate and remove samples that fall 1SD outside of first two principal components   
  vcf.fn <- line
  # Reformat
  snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")
  #snpgdsSummary("test.gds")
  genofile <- snpgdsOpen("test.gds", readonly = F)
  pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)
  pc.percent <- pca$varprop*100
  head(round(pc.percent, 2))
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)
  inliers = tab[tab$EV1 < sd(tab$EV1) & tab$EV1 >-sd(tab$EV1) & tab$EV2 < sd(tab$EV2) & tab$EV2 > -sd(tab$EV2),]
  kep = c(inliers$sample.id)
  closefn.gds(genofile)  

#calculate LD separately for deleterious mutatios (MT=2) and neutral mutations (MT=1)
  multi2 = fread(line)
  df <- multi2 
  df <- lapply(df, gsub, pattern = "1|0", replacement = "1", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "0|1", replacement = "1", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "1|1", replacement = "2", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "0|0", replacement = "0", fixed = TRUE)
  df2 <- data.frame(df)
  
  #deleterious mutations
  dfs = df2[grep("MT=2", df2$INFO),]  
  dfs = dfs[,c(-1:-9)]
  dfs = subset(dfs, select=kep) 
  cols = c(1:(ncol(dfs)))
  dfs[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
  dfs = dfs[(rowSums(dfs)) < 200,]

  
#run this chunk of code to simulate calculatign LD based on polarizing on the reference genome. 
#Ie. if we assume the reference genome is the true ancestral state (but in reality it is most often the major allele, not necessarily the ancestral allele). 
#This chunk of code takes the major alleles in the VCF and specifically sets them as the reference
  dfs$sum = rowSums(dfs)
  norms = dfs[dfs$sum < 100,]
  chans = dfs[dfs$sum > 99,]
  chans <- lapply(chans, gsub, pattern = "0", replacement = "zero", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "2", replacement = "two", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "zero", replacement = "2", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "two", replacement = "0", fixed = TRUE)
  chans <- data.frame(chans)
  cols2s = c(1:(ncol(chans)))
  chans[,cols2s] %<>% lapply(function(x) as.numeric(as.character(x)))
  df3s = rbind(norms,chans)
  df3s$sum = NULL
  df3s$sum = rowSums(df3s)
  df3s = df3s[df3s$sum < 6 ,] #allele coutn cutoff like in the real data, again can be ignored to explore alternative ways of calculating LD
  df3s$sum = NULL
  
  #calculate LD
  mat_s <- data.matrix(df3s) 
  netld_s=(var(colSums(mat_s, na.rm = T)) - sum(rowVars(mat_s, na.rm = T)))/choose(nrow(mat_s), 2)
  burden_s=mean(colSums(mat_s), na.rm = T)
  
  #repeat for neutral mutations
  dfn = df2[grep("MT=1", df2$INFO),]
  dfn = dfn[,c(-1:-9)]
  dfn = subset(dfn, select=kep) 
  cols = c(1:(ncol(dfn)))
  dfn[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
  dfn = dfn[(rowSums(dfn)) < 200,]

  
  dfn$sum = rowSums(dfn)
  normn = dfn[dfn$sum < 100,]
  chann = dfn[dfn$sum > 99,]
  chann <- lapply(chann, gsub, pattern = "0", replacement = "zero", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "2", replacement = "two", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "zero", replacement = "2", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "two", replacement = "0", fixed = TRUE)
  chann <- data.frame(chann)
  cols2n = c(1:(ncol(chann)))
  chann[,cols2n] %<>% lapply(function(x) as.numeric(as.character(x)))
  df3n = rbind(normn,chann)
  df3n$sum = NULL
  df3n$sum = rowSums(df3n)
  df3n = df3n[df3n$sum < 6,]
  df3n$sum = NULL
  
  mat_n <- data.matrix(df3n)  
  netld_n=(var(colSums(mat_n, na.rm = T)) - sum(rowVars(mat_n, na.rm = T)))/choose(nrow(mat_n), 2)
  burden_n=mean(colSums(mat_n), na.rm = T)



  line2=(noquote(paste(netld_s, netld_n, burden_s, burden_n, "1300k_ad"))) #write results to file along with the set if simulations considered (e.g. sims where admixture startes at generation 1300k)
  write(line2,file="../results_100k_netLD.txt",append=TRUE)
}

library(foreach)
library(doParallel)
registerDoParallel(cores=1)
foreach(i=1:nrow(ll), .errorhandling="pass") %dopar% analyze(i)




