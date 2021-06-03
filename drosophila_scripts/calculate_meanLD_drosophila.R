library(data.table)

#very similar to capsella script for calculating mean LD, read in table of genotypes (output from the drosophila_edit_genotypes.sh script)
multi2 = fread("/plas1/george.sandler/dpgp3/ZI_SYN_gts_edited.txt", header = F)

#filter sites with missing data and certain allele frequency cut-offs
multi2$count <- rowSums(is.na(multi2))
multi2 <- multi2[multi2$count < 1,]
multi2$count <-NULL
multi2$count <- rowSums(multi2 ==1, na.rm = T)
multi2 <- multi2[multi2$count < 6,]
multi2 <- multi2[multi2$count > 0,]
multi2$count <-NULL

#convert to matrix and remove non-genotype columns
library(matrixStats)
df2 <- multi2
df2$V1 <- NULL
df2$V2 <- NULL
df2$V3 <- NULL

mat <- data.matrix(df2)

#calulate underdispersion as Sohail et al. 2017 Science
var(colSums(mat, na.rm = T))/sum(rowVars(mat, na.rm = T))

#calculate mean LD per pair of sites
(var(colSums(mat, na.rm = T)) - sum(rowVars(mat, na.rm = T)))/choose(nrow(mat), 2)

