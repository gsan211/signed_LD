#this R script takes as inpu a .txt containing chromosome/site/genotypes columns extracted from a VCF. 
#It shuffles the genotype columns and calculates a NULL distribution for mean LD
#formatted for diploid Capsella dataset

library(matrixStats)
library(data.table)
library(magrittr)
#library(windowscanr)
library(dplyr)

#read in files of interest
multi2 = fread("/plas1/george.sandler/capsella/fly_retry/pop88_LOF_gts.txt")

#remove sites with missing data and sites above a MAC of 5
multi2$count <- rowSums(multi2 == "./.")
multi2 <- multi2[multi2$count <1,]
multi2$count <- NULL
multi2$count <- rowSums(multi2 == "0/1") + 2*rowSums(multi2 == "1/1")
multi2 <- multi2[multi2$count < 6,]
multi2 <- multi2[multi2$count > 0,]
multi2$count <- NULL

#remove no genotype columns and convert genotyeps to derived allelic counts
df <- multi2
df$V1 <- NULL
df$V2 <- NULL
df <- lapply(df, gsub, pattern = "1/0", replacement = "1", fixed = TRUE)
df <- lapply(df, gsub, pattern = "0/1", replacement = "1", fixed = TRUE)
df <- lapply(df, gsub, pattern = "1/1", replacement = "2", fixed = TRUE)
df <- lapply(df, gsub, pattern = "0/0", replacement = "0", fixed = TRUE)
df2 <- data.frame(df)

cols = c(1:(ncol(df2)))
df2[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

#matrix if genotypes
mat <- data.matrix(df2)

#transpose matrix to make shuffling easier
mat2 = t(mat)

#randomly shuffle genotypes WITHIN in each site randomly between individuals
for (i in 1:1000){
  for (i in 1:ncol(mat2)){
    mat2[,i] = sample(mat2[,i], replace=FALSE)
  }

#calculate mean LD (note statistics are reversed compared to the "calculate mean LD" script)
line = ((var(rowSums(mat2, na.rm = T)) - sum(colVars(mat2, na.rm = T)))/choose(ncol(mat2), 2))
write(line,file="/plas1/george.sandler/capsella/fly_retry/review/sig_testing/pop88_LOF_max5_permuted_netLD.txt",append=TRUE)
}

