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

# verification of the paths to the R librairies
.libPaths()

# creation of the directory if it was not created in previous analyses, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create("~/R/x86_64-conda-linux-gnu-library/4.2", recursive = TRUE, showWarnings = FALSE)

# in the case your home directory is not in libPaths
.libPaths('~/R/x86_64-conda-linux-gnu-library/4.2')

# verification of its addition:
.libPaths()

## Code cell 3 ##

# list the required libraries from the CRAN repository
requiredLib <- c(
     "tidyverse",
    "magrittr",
    "matrixStats",
    "BiocManager",
    "ggplot2"
)

# list the required libraries from the Bioconductor project
requiredBiocLib <- c("DESeq2", "vsn")

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

## Code cell 4 ##   

cat("Here is my R session with the loaded packages:\n")
sessionInfo()

## Code cell 5 ##

gohome <- "/shared/projects/2413_rnaseq_cea/"
gohome

myfolder <- getwd()
myfolder

## Code cell 6 ##

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste(myfolder,"/Results/deseq2/", sep = ""), recursive = TRUE)

# storing the path to this output folder in a variable
deseq2folder <- paste(myfolder,"/Results/deseq2/", sep = "")
deseq2folder

# listing the content of the folder
print(system(paste("ls -hlt", deseq2folder), intern = TRUE) )

## Code cell 7 ##

options(repr.plot.width = 15, repr.plot.height = 8) # for figure display in the notebook

## Code cell 8 ##

pca1folder <- paste0(myfolder,"/Results/pca1/")
pca1folder
rdata <- paste0(pca1folder,"RawCounts_Samples.RData")
rdata
load(rdata,verbose = T)
rm(pca1folder) # not used any further

## Code cell 9 ##

ls() 

## Code cell 10 ##

samples # since we have only 11 rows, we can easily print the whole file. We could aslo have used head(samples, n=11).
str(samples)

## Code cell 11 ##

head(countdata)
str(countdata)

## Code cell 12 ##

table(colnames(countdata[,-1]) == samples$SampleName, useNA = "ifany")

## Code cell 13 ##

gencode <- read_tsv("/shared/projects/2413_rnaseq_cea/alldata/Reference/gencode.vM35.primary_assembly.annotation.gtf.gz",
                                     skip = 5, col_names = F, show_col_types = FALSE)


## Code cell 14 ##

str(gencode)
head(gencode)

## Code cell 15 ##

gencode %<>% filter(X3=="gene")
gencode %<>% dplyr::select(X1,X2, X4,X5,X7,X9) 
#gencode %<>% rename(chr=X1, annotation=X2, start=X4, end=X5, strand=X7)
names(gencode)[1:5] <- c("chr", "annotation", "start", "end", "strand")
gencode %<>% separate(X9, into = c("gene_id", "gene_type", "gene_name", "level", "mgi_id", "havana_gene"), sep = ";")
gencode %<>% mutate(gene_id = str_remove(gene_id, "gene_id \""),
                   gene_type = str_remove(gene_type, " gene_type \""), 
                   gene_name = str_remove(gene_name, " gene_name \""),
                   mgi_id = str_remove(mgi_id, " mgi_id \""),
                   havana_gene = str_remove(havana_gene, " havana_gene \""))
gencode %<>% mutate(gene_id = str_remove(gene_id, "\""),
                   gene_type = str_remove(gene_type, "\""), 
                   gene_name = str_remove(gene_name, "\""),
                   mgi_id = str_remove(mgi_id, "\""), 
                   havana_gene = str_remove(havana_gene, "\""))

## Code cell 16 ##

str(gencode)
head(gencode)

## Code cell 17 ##

table(table(gencode$gene_name))
length(unique(gencode$gene_id))
length(unique(gencode$gene_name))

## Code cell 18 ##

write.table(gencode, file = paste0(myfolder, "/Data/gene_info.txt"), quote=F, sep="\t", row.names = F, col.names = T)

## Code cell 19 ##

nrow(countdata)    # 24,432 in counts
str(countdata)

## Code cell 20 ##

row.names(countdata) <- countdata[,"Geneid"] # here we used the name of the column, we could also have put the index (i.e [,1]) without quotes
countdata <- countdata[,-1] # to remove the first column

## Code cell 21 ##

head(countdata)

## Code cell 22 ##

all(names(countdata) == samples$SampleName)

## Code cell 23 ##

row.names(samples) <- samples[,"SampleName"]
head(samples, 11)

## Code cell 24 ##

samples$Condition <- factor(samples$Condition, levels = c("WT", "dHet", "dHetRag" ))

## Code cell 25 ##

table(samples$Condition)

## Code cell 26 ##

levels(samples$Condition)
str(samples$Condition)

## Code cell 27 ##

samples$Condition <- relevel(samples$Condition, ref = "WT")

## Code cell 28 ##

table(samples$Condition)

## Code cell 29 ##

head(samples, n = 2)
str(samples)

## Code cell 30 ##

# round() is used as a security, because DESeq2 works only in integers. Read counts should be integers.
# Be careful: countData is the name of a variable for this function, and countdata is the name of our dataframe

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(countdata, digits = 0),
                              colData = samples,
                              design = ~ Condition)

## Code cell 31 ##

nrow(dds) # 22927 in counts

## Code cell 32 ##
dds

## Code cell 33 ##
head(assay(dds)) # same as coutdata

## Code cell 34 ##

dds2 <- DESeq2::DESeq(dds)

## Code cell 35 ##

dds2

## Code cell 36 ##

head(assay(dds2))

## Code cell 37 ##

head(counts(dds2))

## Code cell 38 ##

head(counts(dds2, normalize = TRUE))

## Code cell 39 ##

data.frame(colData(dds2)[, c("Condition", "Tissue")],
           "sizeFactor" = sizeFactors(dds2),
           "library.size" = colSums(assay(dds2)),
           "total number after norm" = colSums(counts(dds2, normalize = TRUE)))

## Code cell 40 ##

plotDispEsts(dds2)

## Code cell 41 ##

summary(counts(dds2)) # raw counts

## Code cell 42 ##

summary(counts(dds2, normalize = TRUE))# normalized counts

## Code cell 43 ##

# make a colour vector
#conditionColor <- match(samples$Condition, c("dHet", "dHetRag", "WT")) + 1
# '+1' to avoid color '1' i.e. black

# to display the two plots side by side
layout(matrix(1:2, ncol=2))

# Check distributions of samples using boxplots for the non-normalised data in log2
boxplot(log2(counts(dds2)+1), outline = FALSE,
        xlab = "",
        ylab = "Log2(dds2 raw counts)",
        las = 2,
        col = conditionColor,
        main = "Log2(dds2 raw counts)")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(log2(counts(dds2)+1)), col = "blue")

# Check distributions of samples using boxplots
boxplot(log2(counts(dds2, normalize = TRUE)+1), outline = FALSE,
        xlab = "",
        ylab = "DESeq2 norm counts",
        las = 2,
        col = conditionColor,
        main = "DESeq2 norm counts")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median(log2(counts(dds2, normalize = TRUE) + 1)), col = "blue")

## Code cell 44 ##

rlog.dds2 <- rlog(dds2, blind = FALSE)
rlog.dds2

## Code cell 45 ##

# Get log2 counts
counts.rlog.dds2 <- assay(rlog.dds2)

summary(counts.rlog.dds2) # summary for each column

## Code cell 46 ##

# make a colour vector
#conditionColor <- match(samples$Condition, c("dHet", "dHetRag", "WT")) + 1
# '+1' to avoid color '1' i.e. black

# to display the two plots side by side
layout(matrix(1:2, ncol=2))

# Check distributions of samples using boxplots for the non-normalised data in log2
boxplot(log2(counts(dds2)+1), outline = FALSE,
        xlab = "",
        ylab = "Log2(dds2 raw counts)",
        las = 2,
        col = conditionColor,
        main = "Log2(dds2 raw counts)")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(log2(counts(dds2)+1)), col = "blue")

# Check distributions of samples using boxplots for the rlog normalised data
boxplot(counts.rlog.dds2,
        xlab = "",
        ylab = "rlog.dds2 counts",
        las = 2,
        col = conditionColor,
        main = "rlog.dds2 counts")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(counts.rlog.dds2), col = "blue")

## Code cell 47 ##

# Raw counts mean expression Vs standard Deviation (SD)
plot(rowMeans(assay(rlog.dds2)), matrixStats::rowSds(as.matrix(assay(rlog.dds2)), na.rm = TRUE), 
     main = 'rlog normalized read counts: sd vs mean', 
     xlim = c(-1,17),
     ylim = c(-1,5))

## Code cell 48 ##

vsn::meanSdPlot(assay(rlog.dds2), ranks = FALSE)

## Code cell 49 ##

vst.dds2 <- vst(dds2, blind = FALSE)
vst.dds2

## Code cell 50 ##

# Get log2 counts
counts.vst.dds2 <- assay(vst.dds2)

summary(counts.vst.dds2) # summary for each column

## Code cell 51 ##

# make a colour vector
#conditionColor <- match(samples$Condition, c("dHet", "dHetRag", "WT")) + 1
# '+1' to avoid color '1' i.e. black

# to display the two plots side by side
layout(matrix(1:3, ncol = 3))


# Check distributions of samples using boxplots for the non-normalised data in log2
boxplot(log2(counts(dds2)+1), outline = FALSE,
        xlab = "",
        ylab = "Log2(dds2 raw counts)",
        las = 2,
        col = conditionColor,
        main = "Log2(dds2 raw counts)")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(log2(counts(dds2) + 1)), col = "blue")

# Check distributions of samples using boxplots for the rlog normalised data
boxplot(counts.rlog.dds2,
        xlab = "",
        ylab = "rlog.dds2 counts",
        las = 2,
        col = conditionColor,
        main = "rlog.dds2 counts")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(counts.rlog.dds2), col = "blue")

# Check distributions of samples using boxplots
boxplot(counts.vst.dds2,
        xlab = "",
        ylab = "vst.dds2 counts",
        las = 2,
        ylim=c(3,20),
        col = conditionColor,
        main = "vst.dds2 counts")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median.default(counts.vst.dds2), col = "blue")

## Code cell 52 ##

# Raw counts mean expression Vs Standard Deviation (SD)
plot(rowMeans(assay(vst.dds2)), matrixStats::rowSds(as.matrix(assay(vst.dds2)), na.rm = TRUE), 
     main = 'vst normalized read counts: sd vs mean', 
     xlim = c(5,17),
     ylim = c(-1,5))

## Code cell 53 ##

vsn::meanSdPlot(assay(vst.dds2), ranks = FALSE)


## Code cell 54 ##

df <- bind_rows(
  as_data_frame(log2(counts(dds2, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rlog.dds2)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vst.dds2)[, 1:2]) %>% mutate(transformation = "vst"))
  
colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "rlog", "vst")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

## Code cell 55 ##

head(assay(rlog.dds2))

## Code cell 56 ##

rlog.dds2.annot <- merge(as.data.frame(assay(rlog.dds2)),
                gencode,
                by.x = 0, by.y = "gene_id",
                all.x = T, sort = F)
names(rlog.dds2.annot)[1] <- "ensemblID"

## Code cell 57 ##

head(rlog.dds2.annot)

## Code cell 58 ##

resultsNames(dds2)

## Code cell 59 ##

dHet_dHetRag <- c("Condition", "dHet", "dHetRag")

## Code cell 60 ##

res_dHet_dHetRag <- results(dds2, contrast = dHet_dHetRag)
summary(res_dHet_dHetRag)

## Code cell 69 ##

resTot_dHet_dHetRag <- results(dds2, contrast = dHet_dHetRag,  alpha = 0.99999999999999)
summary(resTot_dHet_dHetRag)

## Code cell 70 ##

hist(resTot_dHet_dHetRag$padj, xlim = c(0,1),breaks = 100)

## Code cell 72 ##
hist(resTot_dHet_dHetRag$pvalue, xlim = c(0,0.15),breaks = 200)
hist(resTot_dHet_dHetRag$padj, xlim = c(0,0.15),breaks = 200)

## Code cell 61 ##

res2_dHet_dHetRag <- results(dds2, contrast = dHet_dHetRag,  alpha = 0.05)
summary(res2_dHet_dHetRag)

## Code cell 62 ##

sum(res2_dHet_dHetRag$padj < 0.05, na.rm = TRUE)
summary(res2_dHet_dHetRag)
# 7374

## Code cell 63 ##

res2_dHet_dHetRag_sig_genes <- subset(res2_dHet_dHetRag, padj < 0.05)
dim(res2_dHet_dHetRag_sig_genes)

## 7374     6

## Code cell 64 ##

res2_dHet_dHetRag_sig_ranked <- res2_dHet_dHetRag_sig_genes[order(res2_dHet_dHetRag_sig_genes$pvalue),]
dim(res2_dHet_dHetRag_sig_ranked)
head(res2_dHet_dHetRag_sig_ranked)
## [1] 7374    6

## Code cell 65 ##

res2_dHet_dHetRag_sig_ranked_annot <- merge(as.data.frame(res2_dHet_dHetRag_sig_ranked), gencode, 
                                            by.x = 0, by.y = "gene_id",
                                            all.x = T, sort = F)

## Code cell 66 ##

str(res2_dHet_dHetRag_sig_ranked_annot)

## Code cell 67 ##

names(res2_dHet_dHetRag_sig_ranked_annot)[1] <- "ensemblID"
head(res2_dHet_dHetRag_sig_ranked_annot, n = 20)

## Code cell 68 ##

subset(res2_dHet_dHetRag_sig_ranked_annot, gene_name == "Sos1")

## Code cell 69 ##

write.table(as.data.frame(res2_dHet_dHetRag_sig_ranked_annot),
            file = paste0(deseq2folder,"DESeq2_significant_genes-0_05.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



## Code cell 73 ##

res_dHet_dHetRag_annot <- merge(as.data.frame(res_dHet_dHetRag), gencode, 
                                            by.x = 0, by.y = "gene_id",
                                            all.x = T, sort = F)
str(res_dHet_dHetRag_annot)
names(res_dHet_dHetRag_annot)[1] <- "ensemblID"
head(res_dHet_dHetRag_annot, n = 15)

## Code cell 74 ##

alpha <- 0.00001

## Code cell 75 ##

## Number of genes with an adjusted p-value < alpha *and a positive log(FC)*
dim(res2_dHet_dHetRag_sig_ranked_annot[which(res2_dHet_dHetRag_sig_ranked_annot$padj < alpha & res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange > 2),])[1]
upGenes <- res2_dHet_dHetRag_sig_ranked_annot[which(res2_dHet_dHetRag_sig_ranked_annot$padj < alpha & res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange > 2),]

## Number of genes with an adjusted p-value < alpha *and a negative log(FC)*
dim(res2_dHet_dHetRag_sig_ranked_annot[which(res2_dHet_dHetRag_sig_ranked_annot$padj < alpha & res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange < -2),])[1]
downGenes <- res2_dHet_dHetRag_sig_ranked_annot[which(res2_dHet_dHetRag_sig_ranked_annot$padj < alpha & res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange < -2),]


## Code cell 76 ##

options(repr.plot.width = 10, repr.plot.height = 10)

#draw the plot
plot(-log10(res2_dHet_dHetRag_sig_ranked_annot$padj) ~ res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange, pch = "", xlim = c(-15,15), ylim = c(0,90),
     xlab = "",ylab = "", bty = "n", xaxt = "n", yaxt = "n"  )
title("Leukemic dHet vs non-leukemic dHetRag", font.main = 1, cex.main = 0.9)
axis(1, at = -15:15, tcl = -0.5,cex = 0.7, labels = F )
mtext(-15:15,side = 1,line = 1,at = -15:15, cex = 0.7)
axis(2, at = 0:90, tcl = -0.2, cex = 0.7, labels = F )
mtext(seq(0,90,10),side = 2,line = 0.5, at = seq(0,90,10), cex = 0.7)


## Color in grey genes
points(-log10(res2_dHet_dHetRag_sig_ranked_annot$padj) ~ res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange, pch = 16, cex = 0.5, col = "grey")


## Override the previous colouring for some genes according to the fact they are DE
points(-log10(upGenes$padj) ~ upGenes$log2FoldChange,pch = 16, cex = 0.5, col = "red")
points(-log10(downGenes$padj) ~ downGenes$log2FoldChange, pch = 16,cex = 0.5, col = "green")
mtext("log2(FC)", side = 1, line = 2, cex = 0.8)
mtext("-log10 (adjusted p-value)", side = 2, line = 1.5, cex = 0.8)

abline(h=-log10(alpha), col="blue")
abline(v=-2, col="green")
abline(v=2, col="red")

## Code cell 77 ##

head(res2_dHet_dHetRag_sig_ranked_annot, n = 3)

## Code cell 78 ##

options(repr.plot.width = 10, repr.plot.height = 10)

#draw the plot
plot(-log10(res2_dHet_dHetRag_sig_ranked_annot$padj) ~ res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange, pch = "", xlim = c(-16,16), ylim = c(0,145),
     xlab = "",ylab = "", bty = "n", xaxt = "n", yaxt = "n"  )
title("Leukemic dHet vs non-leukemic dHetRag", font.main = 1, cex.main = 1.5)
axis(1, at = -16:16, tcl = -0.5,cex = 0.7, labels = F )
mtext(-16:16,side = 1,line = 1,at = -16:16, cex = 0.7)
axis(2, at = 0:145, tcl = -0.2, cex = 0.7, labels = F )
mtext(seq(0,145,10),side = 2,line = 0.5, at = seq(0,145,10), cex = 0.7)


## Color in grey genes
points(-log10(res2_dHet_dHetRag_sig_ranked_annot$padj) ~ res2_dHet_dHetRag_sig_ranked_annot$log2FoldChange, pch = 16, cex = 0.5, col = "grey")


## Override the previous colouring for some genes according to the fact they are DE
points(-log10(upGenes$padj) ~ upGenes$log2FoldChange,pch = 16, cex = 0.8, col = "red")
points(-log10(downGenes$padj) ~ downGenes$log2FoldChange, pch = 16,cex = 0.8, col = "blue")
mtext("log2(FC)", side = 1, line = 2, cex = 1.1)
mtext("-log10 (adjusted p-value)", side = 2, line = 1.5, cex = 1.1)

abline(h=-log10(alpha), col="green")
abline(v=-2, col="blue")
abline(v=2, col="red")

## Code cell 79 ##

dim(upGenes)
dim(downGenes)
head(upGenes)

## Code cell 80 ##

write.table(upGenes$gene_name, file=paste0(deseq2folder,"DESeq2_significant_genes-0_00001-up.tsv"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
write.table(downGenes$gene_name, file=paste0(deseq2folder,"DESeq2_significant_genes-0_00001-down.tsv"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)


## Code cell 81 ##

print(ls())

## Code cell 82 ##

save(rlog.dds2.annot, res2_dHet_dHetRag_sig_ranked_annot, gencode, samples, upGenes, downGenes, file = paste0(deseq2folder, "deseq2.RData"))

## Code cell 83 ##

save(res_dHet_dHetRag_annot, file = paste0(deseq2folder, "deseq2_all_genes.RData"))

## Code cell 84 ##   

myfolder
file.copy("/shared/projects/2413_rnaseq_cea/pipeline/Pipe_10-R-Normcounts-exploratory-analysis-II.ipynb", myfolder)


## Code cell 85 ##   

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE, showWarnings = FALSE)

runfolder <- paste0(myfolder,"/run_notebooks")
       
file.copy(paste0(myfolder, "/Pipe_09-R-DESeq2-normalisation-DE.ipynb"), paste0(runfolder, "/Pipe_09-R-DESeq2-normalisation-DE-run.ipynb"))

