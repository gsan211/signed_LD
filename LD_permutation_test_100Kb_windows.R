#Same as other permutation script but calculates LD in 100Kb windows instead
#again formatted to capsella but works for drosophila if the genotype substitution code is changed

library(matrixStats)
library(data.table)
library(magrittr)
library(windowscanr)
library(dplyr)

#read in file with sites of interest
multi2 = fread("/plas1/george.sandler/capsella/fly_retry/pop88_SYN_gts.txt")


multi2$count <- rowSums(multi2 == "./.")
multi2 <- multi2[multi2$count <1,]
multi2$count <- NULL
multi2$count <- rowSums(multi2 == "0/1") + 2*rowSums(multi2 == "1/1")
multi2 <- multi2[multi2$count < 6,]
multi2 <- multi2[multi2$count > 0,]
multi2$count <- NULL

summ = function(X) {M = (sum(X ,na.rm = T)); M}
df3 <- multi2

df3 <- lapply(df3, gsub, pattern = "1/0", replacement = "1", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "0/1", replacement = "1", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "1/1", replacement = "2", fixed = TRUE)
df3 <- lapply(df3, gsub, pattern = "0/0", replacement = "0", fixed = TRUE)
df3 <- data.frame(df3)



df3[,2:ncol(df3)] %<>% lapply(function(x) as.numeric(as.character(x)))
nm <- df3[,c(3:ncol(df3))]

#separate out columns with chromosome/site position so that shufflign doesnt icnlude them
id = df3[,1:2]
id2 = t(id)
df4 = t(nm)
##########################################
##########################################
for (i in 1:1000){
  for (i in 1:ncol(df4)){
    df4[,i] = sample(df4[,i], replace=FALSE)
  }

#after shuffling merge them back together so the window analysis can proceed (requires the first two columns)
df5 = rbind(id2, df4)
df6 = t(df5)
df6 = data.frame(df6)
df6[,2:ncol(df6)] %<>% lapply(function(x) as.numeric(as.character(x)))

#window analysis as in regular LD calculations script
rol_win <- winScan(x = df6,
                   groups = "V1", 
                   position = "V2", 
                   values = c(names(nm)), 
                   win_size = 100000,
                   win_step = 100000,
                   funs = c("summ"))

winn <- rol_win[, -grep("summ$", colnames(rol_win))]
winn <- winn[,c(-1:-4)]
winn <- as.matrix(winn)
a = rowMaxs(winn)

rol_win <- rol_win[, -grep("_n$", colnames(rol_win))]
rol_win <- rol_win[,c(-1:-4)]
mat_win <- as.matrix(rol_win)

line = (var(colSums(mat_win, na.rm = TRUE)) - sum(rowVars(mat_win, na.rm = T)))/sum(combn(a, m = 2, FUN = prod))
write(line,file="/plas1/george.sandler/capsella/fly_retry/permute_inter/pop88_SYN_max5_permuted_netLD100kb.txt",append=TRUE)
}













