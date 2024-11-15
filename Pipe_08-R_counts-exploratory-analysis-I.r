## Code cell 1 ##

session_parameters <- function(){
    
    jupytersession <- c(system('echo "=== Cell launched on $(date) ==="', intern = TRUE),
                        system('squeue -hu $USER', intern = TRUE))
    
    jobid <- system("squeue -hu $USER | awk '/sys/dash {print $1}'", intern = TRUE)
    jupytersession <- c(jupytersession,
                        "=== Current IFB session size: Medium (5CPU, 21 GB) ===",
                        system(paste("sacct --format=JobID,AllocCPUS,ReqMem,NodeList,Elapsed,State -j", jobid), intern = TRUE))
    print(jupytersession[1:6])
    
    return(invisible(NULL))
}

session_parameters()

## Code cell 2 ##

.libPaths()

## Code cell 3 ##

length(list.files("/shared/ifbstor1/software/miniconda/envs/r-4.2.3/lib/R/library"))  # modifier l'indice selon le repertoire souhaité

## Code cell 4 ##

head(installed.packages()[,c(1,2,3)])

grep("ggplot2", installed.packages()[,c(1,2,3)], value = TRUE)

## Code cell 5 ##

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix (thus only created if not yet present)
dir.create("~/R/x86_64-conda-linux-gnu-library/4.2", recursive = TRUE, showWarnings = FALSE)

# this new directory is added to the possible paths:
.libPaths('~/R/x86_64-conda-linux-gnu-library/4.2')

# and we verify its addition:
.libPaths()

## Code cell 6 ##

# list the required libraries from the CRAN repository
requiredLib <- c(
    "tidyverse",
    "data.table",
    "ggfortify",
    "ggrepel",
    "RColorBrewer",
    "matrixStats",
    "BiocManager"    
)

# list the required libraries from the Bioconductor project
requiredBiocLib <- c("affy")

# install required libraries if not yet installed
for (lib in requiredLib) {
  if (!require(lib, character.only = TRUE, quiet = TRUE)) {
    install.packages(lib, quiet = TRUE)
  }
}

for( lib in requiredBiocLib) {
  if (!require(lib, character.only = TRUE, quiet = TRUE)) {
  BiocManager::install(lib, quiet = TRUE)
  }
}

# load libraries
message("Loading required libraries")
for (lib in requiredLib) {
  library(lib, character.only = TRUE)}
for (lib in requiredBiocLib) {
  library(lib, character.only = TRUE)}

# remove variables from the R session if they are no longer necessary 
rm(lib, requiredLib, requiredBiocLib)

## Code cell 7 ##   

cat("---Here are the info of my R session with the loaded packages:\n")
sessionInfo()

## Code cell 8 ##   

getwd()

## Code cell 9 ##

gohome <- "/shared/projects/2413_rnaseq_cea/"
gohome

myfolder <- getwd()
myfolder

## Code cell 10 ##

datafolder <- paste0(myfolder, "/Data/")
datafolder

reffolder <- paste0(gohome, "alldata/Reference/")
reffolder
annot_version <- "vM35"

## Code cell 11 ##

salmonfolder <- paste0(myfolder, "/Results/salmon/")
salmonfolder

## Code cell 12 ##

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/Results/pca1/"), recursive = TRUE)

# storing the path to this output folder in a variable
pca1folder <- paste(myfolder,"/Results/pca1/", sep = "")
pca1folder

# listing the content of the folder
print(system(paste("ls -hlt", pca1folder), intern = TRUE) )

## Code cell 13 ##

options(repr.plot.width = 15, repr.plot.height = 8)

## Code cell 14 ##

file.copy(from = paste0(gohome, "alldata/Data/sampleData-GSE158661.tsv"),
         to = datafolder)

## Code cell 15 ##

samples <- read.table(paste0(datafolder, "sampleData-GSE158661.tsv"),
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE)

## Code cell 16 ##

class(samples)

## Code cell 17 ##

str(samples)

## Code cell 18 ##

head(samples, n = 12)

## Code cell 19 ##

file.copy(from = paste0(gohome, "alldata/Results/featurecounts/11samples_paired-unstranded.counts"),
          to = datafolder)

## Code cell 20 ##

countdata <- read.table(paste0(datafolder, "11samples_paired-unstranded.counts"),
                     sep = "\t",
                     header = TRUE,
                     comment="#")

## Code cell 21 ##

class(countdata)

## Code cell 22 ##

str(countdata)

## Code cell 23 ##

head(countdata)

## Code cell 24 ##

countdata <- countdata[, c(1,7:17)]
str(countdata)
head(countdata)

## Code cell 25 ##

dim(countdata)

## Code cell 26 ##

names(countdata)

## Code cell 27 ##

# We remove the common prefix in all columns except the first one
colnames(countdata)[-1] <- gsub("X.shared.projects.2413_rnaseq_cea.alldata.Results.featurecounts.",
                         "",
                         colnames(countdata)[-1])

# We remove the common suffix in all columns except the first one
colnames(countdata)[-1] <- gsub("_Aligned.sortedByNames.bam",
                         "",
                         colnames(countdata)[-1])

## Code cell 28 ##

names(countdata)

## Code cell 29 ##

table(colnames(countdata)[-1] == samples$SampleID)

## Code cell 30 ##

colnames(countdata)[-1] <- samples$SampleName
str(countdata)
head(countdata)

## Code cell 31 ##

cat("Start time is :", date(), "\n")
cat("Gencode reference : ", annot_version, "\n")  
cat("Salmon directory is ", salmonfolder, "\n")

## Code cell 32 ##

R.utils::copyDirectory(from = paste0(gohome, "alldata/Results/salmon/"),
                       to = salmonfolder,
                       recursive = TRUE)

## Code cell 33 ##

sample_names <- list.files("./Results/salmon/")
sample_names

## Code cell 34 ##

sample_files <- paste0("./Results/salmon/",sample_names, "/quant.sf")
sample_files

## Code cell 34 ##

cat("\nNumber of samples : ", length(sample_names), "\n")

## Code cell 36 ##

cat("\nFirst sample :", sample_files[1], "\n")

## Code cell 37 ##

numreads <- read.table(sample_files[1], sep = "\t", header= T)
numreads <- numreads [,c(1,5)]
head(numreads)
dim(numreads)

## Code cell 38 ##

for (i in 2:length(sample_files)) { 
 mytab <- read.table(sample_files[i], sep = "\t", header= T)
 mytab  <- mytab  [,c(1,5)]  
 names(mytab)[2] <- sample_names[i]
 numreads <- numreads %>% left_join(mytab, by='Name')
}
names(numreads) <- c("transcript_id", sample_names)	
head(numreads)
rm(mytab)

## Code cell 39 ##

table(names(numreads)[-1] == samples$SampleID)
names(numreads)[-1]  <- samples$SampleName
head(numreads)

## Code cell 40 ##

tpm <- read.table(sample_files[1], sep = "\t", header= T)
tpm <- tpm [,c(1,4)]
#head(tpm)

for (i in 2:length(sample_files)) { 
 mytab <- read.table(sample_files[i], sep = "\t", header= T)
 mytab  <- mytab  [,c(1,4)]  
 names(mytab)[2] <- sample_names[i]
 tpm <- tpm %>% left_join(mytab, by='Name')
}
names(tpm) <- c("transcript_id", sample_names)	
cat("Table dimensions : ", dim(tpm), "\n")
head(tpm)
rm(mytab)

## Code cell 41 ##

table(names(tpm)[-1] == samples$SampleID)
names(tpm)[-1]  <- samples$SampleName
head(tpm)

## Code cell 42 ##

transcripts_info <- data.table::fread(paste0(reffolder, "salmon/gencode.vM35.transcripts.fa.gz"), header = FALSE) %>% filter(str_detect(V1, "^>"))
str(transcripts_info)

## Code cell 43 ##

transcripts_info <- transcripts_info %>% separate(V1,
                                                  into = c("transcript_id", "gene_id", "havana_gene_id", "havana_tr_id", "ext_tr_name", "gene_name", "length", "biotype"),
                                                  remove = TRUE,
                                                 sep = "\\|",
                                                 extra ="drop")
str(transcripts_info)

## Code cell 44 ##

transcripts_info <- transcripts_info %>% mutate(transcript_id = str_replace(transcript_id, ">", ""))
head(transcripts_info)
tail(transcripts_info)

## Code cell 45 ##

cat("\nCheck missing transcripts...\n")
cat(" this experiment - the reference :")
cat(length(setdiff(numreads$transcript_id, transcripts_info$transcript_id)), "\n")
# [1] 0
cat(" the reference - this experiment :")
cat(length(setdiff(transcripts_info$transcript_id, numreads$transcript_id)), "\n")
# [1] 847


## Code cell 46 ##

genecounts <- numreads  %>% left_join(transcripts_info %>% select(transcript_id, gene_id),
                                      by='transcript_id') %>% 
  select(transcript_id, gene_id, everything()) %>% select(-1) %>% 
  group_by(gene_id) %>% summarize_all(sum, na.rm = TRUE) 

#genecounts <- genecounts %>% mutate(gene_id = str_replace(gene_id, "\\.\\d{1,2}$", "")) #in case you want to remove the .XX in the gene_id

cat("Table dimensions:", dim(genecounts), "\n")

head(genecounts)

## Code cell 47 ##

genecounts_symbols <- genecounts %>% left_join(transcripts_info %>%
                                   mutate(gene_id = str_remove(gene_id, ".\\d{1,2}$")) %>%  
 select(gene_id, gene_name), by = 'gene_id') %>% select(-gene_id) %>% 
        select(gene_name, everything()) %>% distinct()

cat("Table dimensions:", dim(genecounts_symbols), "\n")

head(genecounts_symbols)

## Code cell 48 ##

gtpm <- tpm  %>% left_join(transcripts_info %>% select(transcript_id, gene_id), by='transcript_id') %>% 
	select(transcript_id, gene_id, everything()) %>% select(-1) %>% 
	group_by(gene_id) %>% summarize_all(sum, na.rm = TRUE) 

#gtpm <- gtpm %>% mutate(gene_id = str_replace(gene_id, "\\.\\d{1,2}$", "")) #in cas eyou want to remove the .XX in the gene_id

cat("Table dimensions:", dim(gtpm), "\n")

head(gtpm)

## Code cell 49 ##

gtpm_symbols <- gtpm %>% left_join(transcripts_info %>%
                                   mutate(gene_id = str_remove(gene_id, ".\\d{1,2}$")) %>%  
 select(gene_id, gene_name), by = 'gene_id') %>% select(-gene_id) %>% 
        select(gene_name, everything()) %>% distinct()

cat("Table dimensions:", dim(gtpm_symbols), "\n")

head(gtpm_symbols)

## Code cell 50 ##

cat("End time is :", date(), "\n")

## Code cell 51 ##

# keeping outcome in a vector of 'logical' values (ie TRUE or FALSE, or NA)
keep <- rowSums(countdata[-1]) > 10

# summary of test outcome: number of genes in each class:
table(keep, useNA="always") 

## Code cell 52 ##


# subset genes where test was TRUE
countdata <- countdata[keep,]
rm(keep)

# check dimension of new count matrix
dim(countdata)

# 24432 12

## Code cell 53 ##

# keeping outcome in a vector of 'logical' values (ie TRUE or FALSE, or NA)
keep <- rowSums(genecounts[-1]) > 10

# summary of test outcome: number of genes in each class:
table(keep, useNA="always") 

## Code cell 54 ##


# subset genes where test was TRUE
genecounts <- genecounts[keep,]
rm(keep)

# check dimension of new count matrix
dim(genecounts)

# 21427 12

## Code cell 55 ##

cat("----number of common genes:\n")

shared_star_salmon <- intersect(genecounts$gene_id, countdata$Geneid)
length(shared_star_salmon)

cat("----number of genes only in STAR:\n")
length(setdiff( countdata$Geneid, genecounts$gene_id))

cat("----number of genes only in Salmon:\n")
length(setdiff(genecounts$gene_id, countdata$Geneid))

## Code cell 56 ##

summary(countdata[,-1])

## Code cell 57 ##

summary(genecounts[,-1])

## Code cell 58 ##

opar <- par(no.readonly=TRUE) # l'argument no.readonly=TRUE permet de supprimer l'affichage d'éventuels warnings 
par(mar = c(8,2,2,2)) # to increase margin at the bottom to display full sample names
boxplot(countdata[,-1], main = 'STAR/featurecounts raw readcounts distribution across samples', xaxt = "n") # no display of xlabels
text(1:11, y = par("usr")[3] - 0.45,
     labels = names(countdata)[-1],
     srt = 30, adj = 1, xpd = NA, cex= 0.9)# add xlabels with 30 degres srotation
suppressWarnings(par(opar))


## Code cell 59 ##

opar <- par(no.readonly=TRUE) # l'argument no.readonly=TRUE permet de supprimer l'affichage d'éventuels warnings 
par(mar = c(8,2,2,2)) # to increase margin at the bottom to display full sample names
boxplot(genecounts[,-1], main = 'Salmon gene raw readcounts distribution across samples', xaxt = "n") # no display of xlabels
text(1:11, y = par("usr")[3] - 0.45,
     labels = names(countdata)[-1],
     srt = 30, adj = 1, xpd = NA, cex= 0.9)# add xlabels with 30 degres srotation
suppressWarnings(par(opar))


## Code cell 60 ##

opar <- par(no.readonly=TRUE)
par(mfrow = c(1,2))
# Raw counts mean expression vs Standard Deviation (SD) with STAR/featurecounts
plot(rowMeans(countdata[,-1]), apply(countdata[,2:12], 1, sd, na.rm = TRUE), 
     main = 'STAR/feature counts Raw read counts: sd vs mean', 
     xlim = c(0,100000), # to zoom
     ylim = c(0,50000), # to zoom
     xlab = "Mean of raw counts per gene",
     ylab = "SD of raw counts per gene"
    )
# Raw counts mean expression vs Standard Deviation (SD) with Salmon
plot(rowMeans(genecounts[,-1]), apply(genecounts[,2:12], 1, sd, na.rm = TRUE), 
     main = 'Salmon raw read counts: sd vs mean', 
     xlim = c(0,100000), # to zoom
     ylim = c(0,50000), # to zoom
     xlab = "Mean of raw counts per gene",
     ylab = "SD of raw counts per gene"
    )
par(opar) 

## Code cell 61 ##

summary(log2(countdata[,2]+1)) # summary for first sample column 2

## Code cell 62 ##

summary(log2(countdata[,2:12]+1)) # summary for each sample

## Code cell 63 ##

summary(log2(genecounts[,2:12]+1)) # summary for each sample

## Code cell 64 ##

# make a colour vector
conditionColor <- match(samples$Condition, c("dHet", "dHetRag", "WT")) + 1
# '+1' to avoid color '1' i.e. black

# Check distributions of samples using boxplots
opar <- par(no.readonly=TRUE) # l'argument no.readonly=TRUE permet de supprimer l'affichage d'éventuels warnings 
par(mar = c(6,2,2,2))
boxplot(log2(countdata[,2:12]+1),
        xlab = "", xaxt = "n",
        ylab = "Log2(Counts)",
        las = 2,
        col = conditionColor,
        main = "Log2( STAR/featurecounts Counts) distribution across samples")
text(1:11, y = par("usr")[3] - 0.45,
     labels = names(countdata)[-1],
     srt = 30, adj = 1, xpd = NA, cex= 0.9)


# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(as.matrix(log2(countdata[,2:12]+1))), col="blue")

suppressWarnings(par(opar))

## Code cell 65 ##

# Check distributions of samples using boxplots
opar <- par(no.readonly=TRUE) # l'argument no.readonly=TRUE permet de supprimer l'affichage d'éventuels warnings 
par(mar = c(6,2,2,2))
boxplot(log2(genecounts[,2:12]+1),
        xlab = "", xaxt = "n",
        ylab = "Log2(Counts)",
        las = 2,
        col = conditionColor,
        main = "Log2(STAR Counts) distribution across samples")
text(1:11, y = par("usr")[3] - 0.45,
     labels = names(genecounts)[-1],
     srt = 30, adj = 1, xpd = NA, cex= 0.9)

# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(as.matrix(log2(genecounts[,2:12]+1))), col="blue")

suppressWarnings(par(opar))

## Code cell 66 ##

# Log2 counts standard deviation (sd) vs mean expression

opar <- par(no.readonly=TRUE)
par(mfrow = c(1,2))
plot(rowMeans(log2(countdata[,2:12]+1)),
              matrixStats::rowSds(as.matrix(log2(countdata[,2:12]+1), na.rm = TRUE)), 
     main = 'Log2 STAR/Featurecounts Counts: sd vs mean',
     xlab = "Mean of log2(raw counts) per gene",
     ylab = "SD of log2(raw counts) per gene"
    
    )
plot(rowMeans(log2(genecounts[,2:12]+1)),
              matrixStats::rowSds(as.matrix(log2(genecounts[,2:12]+1), na.rm = TRUE)), 
     main = 'Log2 STAR/Featurecounts Counts: sd vs mean',
     xlab = "Mean of log2(raw counts) per gene",
     ylab = "SD of log2(raw counts) per gene"
    
    )
par(opar) 

## Code cell 67 ##

tlogcounts_star <- t(log2(countdata[,2:12]+1))
dim(tlogcounts_star)

## Code cell 68 ##

# run PCA
PCAdata <- prcomp(tlogcounts_star)

# plot PCA
autoplot(PCAdata)


## Code cell 69 ##

autoplot(PCAdata,
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
    geom_text_repel(aes(x = PC1, y = PC2, label = SampleName),
                        box.padding = 0.8)


## Code cell 70 ##

autoplot(PCAdata,
         x = 2,    # PC2
         y = 3,    # PC3
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
    geom_text_repel(aes(x = PC2, y = PC3, label = SampleName),
                    box.padding = 0.8)


## Code cell 71 ##

autoplot(PCAdata,
         x = 3,    # PC3
         y = 4,    # PC4
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
    geom_text_repel(aes(x = PC3, y = PC4, label = SampleName),
                    box.padding = 0.8)

rm(PCAdata) # we remove the PCA dtaa from the session for memory reasons

## Code cell 72 ##

tlogcounts_salmon <- t(log2(genecounts[,2:12]+1))

PCAdata <- prcomp(tlogcounts_star)

autoplot(PCAdata,
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
    geom_text_repel(aes(x = PC1, y = PC2, label = SampleName),
                        box.padding = 0.8)


## Code cell 73 ##

clusters <- hclust(dist(as.matrix(tlogcounts_star)), method ="ward.D")
plot(clusters, labels = samples$SampleName)

rm(clusters, tlogcounts_star)    

## Code cell 74 ##

clusters <- hclust(dist(as.matrix(tlogcounts_salmon)), method ="ward.D")
plot(clusters, labels = samples$SampleName)

rm(clusters, tlogcounts_salmon)    

## Code cell 75 ##

affy::plotDensity(log2(countdata[,2:12]+1),
                  xlab("Density"),ylab("log2(Counts)"),
                  col = 1:11)
legend(x = 14, y = 0.13,legend = names(log2(countdata[,2:12]+1)),
       col = 1:11, lty = 1:11, bty = "n")

## Code cell 76 ##

affy::plotDensity(log2(genecounts[,2:12]+1),
                  xlab("Density"),ylab("log2(Counts)"),
                  col = 1:11)
legend(x = 14, y = 0.13,legend = names(log2(genecounts[,2:12]+1)),
       col = 1:11, lty = 1:11, bty = "n")

## Code cell 77 ##

temp1 <- pivot_longer(subset(countdata, Geneid %in% shared_star_salmon),cols=2:12, names_to = "sample", values_to = "counts")
temp2 <- pivot_longer(subset(genecounts, gene_id %in% shared_star_salmon),cols=2:12, names_to = "sample", values_to = "counts")
temp1 <- temp1[order(temp1$Geneid),]
temp2 <- temp2[order(temp2$gene_id),]
cor.test(temp1$counts, temp2$counts)
plot(log2(temp1$counts+1), log2(temp2$counts+1))

rm(temp1, temp2)

## Code cell 78 ##  

ls()
save(countdata, genecounts, genecounts_symbols, gtpm, gtpm_symbols, samples, conditionColor, file = paste0(pca1folder,"RawCounts_Samples.RData"))
## Code cell 79 ##   

myfolder
file.copy("/shared/projects/2413_rnaseq_cea/pipeline/Pipe_09-R-DESeq2-normalisation-DE.ipynb", myfolder)

## Code cell 80 ##   

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE)

runfolder <- paste0(myfolder,"/run_notebooks")
       
# file.copy(paste0(myfolder, "/Pipe_08-R_counts-exploratory-analysis-I.ipynb"), runfolder)
file.copy(paste0(myfolder, "/Pipe_08-R_counts-exploratory-analysis-I.ipynb"), paste0(runfolder, "/Pipe_08-R_counts-exploratory-analysis-I-run.ipynb"))
