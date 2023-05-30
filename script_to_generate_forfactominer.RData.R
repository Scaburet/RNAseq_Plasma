#!user/bin/env Rscript

##########################################################
### Script to prepare the data for the factominer analysis
## claire.vandiedonck@univ-paris-diderot.fr
##########################################################

#------------------
#----- import expression data:
#------------------

setwd(".")# modify the path if necessary
load("norm.quant.Rdata")
ls()
# [1] "norm.quant"  
str(norm.quant)
# num [1:47323, 1:264] 7.35 7.09 6.3 6.57 6.49 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:47323] "6450255" "2570615" "6370619" "2600039" ...
# ..$ : chr [1:264] "5753669129_B" "5753669129_C" "5753669129_E" "5753669129_K" ...
class(norm.quant)
#[1] "matrix"

#------------------
#-----import metadata:
#------------------

load("samples_info.RData")
ls()
# "norm.quant"   "samples_info"
str(samples_info)
# 'data.frame':	264 obs. of  8 variables:
#   $ array.labels: chr  "5753669129_B" "5753669129_C" "5753669129_E" "5753669129_K" ...
# $ PedID       : chr  "1" "1" "1" "1" ...
# $ ID          : chr  "1" "1" "1" "1" ...
# $ Status      : chr  "2" "2" "2" "2" ...
# $ Stim        : chr  "0" "0" "6" "6" ...
# $ Full        : chr  "1_1_2_0" "1_1_2_0" "1_1_2_6" "1_1_2_6" ...
# $ Sex         : chr  "1" "1" "1" "1" ...
# $ Age         : int  10 10 10 10 10 10 18 18 18 18 ...

#------------------
#----- get gene list according to their symbol from probes data object:
#------------------

load("probes.Rdata")  
ls()
# [1] "norm.quant"   "probes"       "samples_info"
str(probes)
# 'data.frame':	47323 obs. of  9 variables:
# $ ProbeID              : int  6450255 2570615 6370619 2600039 2650615 5340672 2000519 3870044 7050209 1580181 ...
# $ CHROMOSOME           : chr  "7" "19" "19" "10" ...
# $ CYTOBAND             : chr  "7p15.3e" "19q13.43c" "19q13.43c" "10q11.23c" ...
# $ PROBE_CHR_ORIENTATION: chr  "-" "-" "-" "-" ...
# $ PROBE_COORDINATES    : chr  "20147187-20147236" "63548541-63548590" "63549180-63549229" "52566586-52566635" ...
# $ REFSEQ_ID            : chr  "NM_182762.2" "NM_130786.2" "NM_130786.2" "NM_138932.1" ...
# $ ENTREZ_GENE_ID       : int  346389 1 1 29974 29974 29974 23784 23784 23784 54715 ...
# $ TargetID             : chr  "7A5" "A1BG" "A1BG" "A1CF" ...
# $ SYMBOL               : chr  "7A5" "A1BG" "A1BG" "A1CF" ...

gene_list <- probes$SYMBOL
length(gene_list)
# [1] 47323

length(which(is.na(gene_list)))
#[1] 0
length(gene_list[gene_list==""])
#[1] 3270 # there is no matching even in biomaRt for the corresponding probes -> we keep the probeID in the column
gene_list[gene_list==""] <- paste0("probe_", probes$ProbeID[which(probes$SYMBOL=="")])
head(gene_list)
# [1] 7A5  A1BG A1BG A1CF A1CF A1CF

#------------------
#-----compute the mean of the probes for the same gene and put it in a matrix:
#------------------

gene_average <- tapply(norm.quant[,1], gene_list, mean, na.rm=TRUE) # le nom de la ligne est cree automatiquement par tapply()
for (i in 2:ncol(norm.quant)){
  this_gene <- tapply(norm.quant[,i], gene_list, mean, na.rm=TRUE)
  gene_average <- cbind(gene_average, this_gene)
}
dim(gene_average)
# [1] 34696   264
class(gene_average)
# [1] "matrix"

# get sample ids and name columns gene_average
colnames(gene_average) <- samples_info$Full

#------------------
#-----compute the mean of biological relicates
#------------------

# not all samples are duplicated, there are some singletons:
# extract singletons:
singletons <- names(which(table(colnames(gene_average))==1))
singletons
# [1] 2_5_1_0  4_10_1_0 5_13_1_24  7_16_1_0  8_20_1_0  8_21_1_6 

gene_average_dup <- gene_average[,!colnames(gene_average) %in% singletons]
dim(gene_average_dup)
#[1] 34696   258 # the 6 singletons are no longer present

# compute mean on duplicates:
duplicates_means <- apply(gene_average_dup[,1:2], 1, mean, na.rm=T)
length(duplicates_means)
# [1] 34696
for (i in seq(3,ncol(gene_average_dup),2)){
  this_dup <- apply(gene_average_dup[,i:(i+1)], 1, mean)
  duplicates_means <- cbind(duplicates_means, this_dup)
}
dim(duplicates_means) 
# [1] 34696   129
class(duplicates_means)
# [1] "matrix"
colnames(duplicates_means) <- unique(colnames(gene_average_dup))

# we add singletons in the same object:
duplicates_means <- cbind(duplicates_means, gene_average[,singletons])
dim(duplicates_means)
#[1] 34696   135

# to finish reorder columns as in the original dataframe
duplicates_means <- duplicates_means[,unique(samples_info$Full)]
class(duplicates_means)
#[1] "matrix"

# duplicates_means is the matrix of the data after averaging acroos probes and biological replicates
save(duplicates_means, file="duplicate_means.RData")

#------------------
#-----the most variable genes:
#------------------

#We will use 1920 obestrevations for Factominer

# write a generic function to compute the coefficient of variation, which is the sdandard deviation divided by the mean:
compute_variation_coef <- function(x){
  coef_var <- sd(x, na.rm=T)/mean(x, na.rm=T)
  return(coef_var)
}

# apply the above function to the gene mean expression levels
mycv <- apply(duplicates_means, 1, compute_variation_coef)

#check the distribution of the coefficients of variation:
summary(mycv)
quantile(mycv, seq(0,1, 0.1))
# 0%         10%         20%         30%         40%         50%         60%         70%         80% 
# 0.003407795 0.012913757 0.015161366 0.016234326 0.017105226 0.018033082 0.019232827 0.021150973 0.024767016 
# 90%        100% 
#   0.032353025 0.282272056

# identify the 1920 most variable genes
top_cv <- sort(mycv, decreasing = T)[1:1920]
# whoch proportion of genes do we keep?
1920/length(mycv)
#[1] 0.05533779 # it correspond to the top 5.5% variable genes
# what is the minimal coefficient of variation we keep?
top_cv[1920]
# INSL3 
# 0.0406962 # the min coef of variation kept for factominer
  
# make a dataframe subset keeping only the 1920 most variable genes based on their coefficient of variation
HighestCV <- duplicates_means[row.names(duplicates_means) %in% names(top_cv),]
dim(HighestCV)
#[1] 1920  135
class(HighestCV)
# [1] "matrix"
save(HighestCV, file="HighestCV.RData")

#------------------
#-----finalyze the dataframe for factominer analysis:
#------------------

# get the info from the 135 unique samples (duplicates counted once + singletons)
unique_samples_infos <- unique(samples_info[,-1])
row.names(unique_samples_infos) <- unique_samples_infos$Full
unique_samples_infos <- unique_samples_infos[,-5]
str(unique_samples_infos)
# 'data.frame':	135 obs. of  6 variables:
# $ PedID : chr  "1" "1" "1" "1" ...
# $ ID    : chr  "1" "1" "1" "2" ...
# $ Status: chr  "2" "2" "2" "1" ...
# $ Stim  : chr  "0" "6" "24" "0" ...
# $ Sex   : chr  "1" "1" "1" "1" ...
# $ Age   : int  10 10 10 18 18 18 14 14 14 10 ...
  
# generate on dataframe for factominer with the mean values in the first columns, then the sample infos, each sample beaing in one row  
for_factominer <- data.frame(t(HighestCV), unique_samples_infos, stringsAsFactors = F)# si tu veux garder la factorisation, ne mets pas le dernier argument
# le format est un dataframe est non une matrice qui sinon aurait converti toutes les valeurs numeriques en chaines de caracteres
dim(for_factominer)
#[1]  135 1926

save(for_factominer, file="for_factominer.RData")
  