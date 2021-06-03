library(matrixStats)
library(data.table)
library(magrittr)
library(windowscanr)
library(dplyr)

#read in file with site information and genotypes of individuals e.g.
multi2 = fread(../pop88_SYN_gts.txt)

#count number of individuals with no genotypes info, remove all sites whith missing info
multi2$count <- rowSums(multi2 == "./.")
multi2 <- multi2[multi2$count <1,]
multi2$count <- NULL

#count number of derived alleles at each site, filter out sites with more than 5 copies or less than 1
multi2$count <- rowSums(multi2 == "0/1") + 2*rowSums(multi2 == "1/1")
multi2 <- multi2[multi2$count < 6,]
multi2 <- multi2[multi2$count > 0,]
multi2$count <- NULL

#remove columns with non-genotypic data, convert genotypes to #derived alleles (i.e. heterozygotes -->1, homozygotes --> 2 etc.)
df <- multi2
df$V2 <- NULL
df$V1 <- NULL
df <- lapply(df, gsub, pattern = "1/0", replacement = "1", fixed = TRUE)
df <- lapply(df, gsub, pattern = "0/1", replacement = "1", fixed = TRUE)
df <- lapply(df, gsub, pattern = "1/1", replacement = "2", fixed = TRUE)
df <- lapply(df, gsub, pattern = "0/0", replacement = "0", fixed = TRUE)
df2 <- data.frame(df)

#convert to numeric
cols = c(1:(ncol(df2)))
df2[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

#convert to matrix for variance calculations
mat <- data.matrix(df2)

#calculate mean LD per pair of mutations
(var(colSums(mat, na.rm = T)) - sum(rowVars(mat, na.rm = T)))/choose(nrow(mat), 2)

#alternatively calculate underdispersion (like the stats presented in Sohail et al. 2017, Science)
(var(colSums(mat, na.rm = T))/sum(rowVars(mat, na.rm = T)))

####################################################################################################
####################################################################################################
#calculate mean LD in 100Kb windows


summ = function(X) {M = (sum(X ,na.rm = T)); M}
df3 <- multi2

#same as before, replace genotype calls with #derived alleles
df3 <- lapply(df3, gsub, pattern = "1/0", replacement = "1", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "0/1", replacement = "1", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "1/1", replacement = "2", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "0/0", replacement = "0", fixed = TRUE)
df3 <- data.frame(df3)


#convert to numeric and create a dataframe with only the columns containin genotypic information (shortcut for the window analysis coming up)
df3[,2:ncol(df3)] %<>% lapply(function(x) as.numeric(as.character(x)))
nm <- df3[,c(3:ncol(df3))]

#rolling window analysis using the windowscanr package
#V1 = chromosome column, V2 = position, values = all the genotypic columns (defined by the columns present in dataframe nm)

rol_win <- winScan(x = df3,
                   groups = "V1", 
                   position = "V2", 
                   values = c(names(nm)), 
                   win_size = 100000,
                   win_step = 100000,
                   funs = c("summ"))

#removing output columns that aren't of use
winn <- rol_win[, -grep("summ$", colnames(rol_win))]
winn <- winn[,c(-1:-4)]
winn <- as.matrix(winn)

#calculate mean LD in 100Kb windows as descriped in methods of Sandler at al. 2021 MBE
a = rowMaxs(winn)
rol_win <- rol_win[, -grep("_n$", colnames(rol_win))]
rol_win <- rol_win[,c(-1:-4)]
mat_win <- as.matrix(rol_win)

#underdispersion LD statistic
var(colSums(mat_win, na.rm = TRUE))/sum(rowVars(mat_win, na.rm = T))

#mean LD per pair of sites
(var(colSums(mat_win, na.rm = TRUE)) - sum(rowVars(mat_win, na.rm = T)))/sum(combn(a, m = 2, FUN = prod))
