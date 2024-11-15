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


# list the required libraries from the CRAN repository
requiredLib <- c(
    "ggfortify",
    "ggrepel",
    "RColorBrewer",
    "ggplot2",
    "stringr",
    "matrixStats",
    "BiocManager",
    "ggnewscale"
)

# list the required libraries from the Bioconductor project
requiredBiocLib <- c("DESeq2",
                     "org.Mm.eg.db",
                     "clusterProfiler",
                     "enrichplot",
                     "ComplexHeatmap",
                     "ReactomePA")

# install required libraries if not yet installed
for (lib in requiredLib) {
  if (!require(lib, character.only = TRUE, quiet = TRUE)) {
    install.packages(lib, quiet = TRUE, repos = "https://cloud.r-project.org")
  }
}
for( lib in requiredBiocLib) {
  if (!require(lib, character.only = TRUE, quiet = TRUE)) {
  BiocManager::install(lib, quiet = TRUE, update = FALSE)
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

## Code cell 3 ##   

cat("Here is my R session with the loaded packages:\n")
sessionInfo()

## Code cell 4 ##


gohome <- "/shared/projects/2413_rnaseq_cea/"
gohome

myfolder <- getwd()
myfolder


## Code cell 5 ##

# creation of the output directories, recursive = TRUE is equivalent to the mkdir -p in Unix
# and will generate a warning if the directory already exists

dir.create(paste0(myfolder,"/Results/enrich/"), recursive = TRUE)
dir.create(paste0(myfolder,"/Results/gsea/"), recursive = TRUE)
dir.create(paste0(myfolder,"/Results/wgcna/"), recursive = TRUE)

# storing the paths to the output folders in variables
enrichfolder <- paste0(myfolder,"/Results/enrich/")
enrichfolder

gseafolder <- paste0(myfolder,"/Results/gsea/")
gseafolder

wgcnafolder <- paste0(myfolder,"/Results/wgcna2/")
wgcnafolder

# specifying the path to the folder with input data
deseq2folder <- paste0(myfolder,"/Results/deseq2/")
deseq2folder

# listing the content of the folder
print(system(paste("ls -hlt", deseq2folder), intern = TRUE) )

## Code cell 6 ##

options(repr.plot.width = 15, repr.plot.height = 10) # for figure display in the notebook

## Code cell 7 ##

ls()

## Code cell 8 ##

rdata <- paste0(deseq2folder, "deseq2-final.RData")
rdata
load(rdata, verbose = TRUE)

## Code cell 9 ##

deseq2folder
gohome

## Code cell 10 ##

ls()

## Code cell 11 ##

str(top50DE)

## Code cell 12 ##

str(downGenes)

## Code cell 13 ##

str(upGenes)

## Code cell 14 ##

print(list.files("./Results/deseq2/", pattern = "DESeq2_significant_genes-0_00001"))

## Code cell 15 ##

rdata <- paste0(deseq2folder, "deseq2_all_genes.RData")
rdata
load(rdata, verbose = TRUE)

## Code cell 16 ##

rm(conditionColor, dup_genes, eig.val, for_factominer, gencode,
  PCAres, rdata, res2_dHet_dHetRag_sig_ranked_annot, top50var, var_genes)  
ls()

# Code cell 17 ##

mygenes <- unique(top50DE$gene_name) 

# Code cell 18 #

enr_go <- enrichGO(gene = mygenes, 
             ont ="ALL", 
             OrgDb = org.Mm.eg.db,
             universe = unique(res_dHet_dHetRag_annot$gene_name),
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.1, 
             pAdjustMethod = "BH")

# Code cell 19 #

head(enr_go, 10) # but only 7 are significant!

# Code cell n°20

options(repr.plot.width = 12, repr.plot.height = 7)
barplot(enr_go, showCategory = 7,
        label_format = function(x) stringr::str_wrap(x, width = 120)) + # to be able to see terms description in a single row : play on the number (eg. 120)
        ggtitle("barplot for ORA") ## uses ggplot2 

# Code cell n°21

options(repr.plot.width = 12, repr.plot.height = 7)
dotplot(enr_go,
        showCategory = 9)
        label_format = function(x) stringr::str_wrap(x, width=120) +
        ggtitle("dotplot for ORA")

# Code cell n°22

options(repr.plot.width = 30, repr.plot.height= 10)
enr_go <- enrichplot::pairwise_termsim(enr_go, method = "JC", semData = "org.Mm.eg.db")
emapplot(enr_go, showCategory = 20)
# Raw cell 1 ##

options(repr.plot.width = 30, repr.plot.height = 15)
treeplot(enr_go,
         hclust_method = "average")
# Code cell n°23

options(repr.plot.width = 15, repr.plot.height = 10)

upsetplot(enr_go)

# Code cell n°24

cnetplot(enr_go, categorySize ="pvalue", showCategory = 5)

# Code cell n°25

mygenes_entrezid <- clusterProfiler::bitr(mygenes,
                                          fromType = 'SYMBOL',
                                          toType = 'ENTREZID',
                                          OrgDb = 'org.Mm.eg.db')
mygenes_entrezid

# Code cell n°26

enr_react <- enrichPathway(gene = mygenes_entrezid$ENTREZID,
                           organism = "mouse",
                           pvalueCutoff = 0.20,
                           readable = TRUE)
head(enr_react)


# Code cell n°27

options(repr.plot.width = 12, repr.plot.height = 7)
dotplot(enr_react,
        showCategory = 5)+
        ggtitle("dotplot for ORA")

# Code cell n°28

mygenes_uniprot <- clusterProfiler::bitr(mygenes,
                                         fromType = 'SYMBOL',
                                         toType = 'UNIPROT',
                                         OrgDb = 'org.Mm.eg.db')
mygenes_uniprot

# Code cell n°29

table(mygenes_uniprot$SYMBOL)

# Code cell n°30

mygenes_uniprot <- mygenes_uniprot[!duplicated(mygenes_uniprot$SYMBOL), ]
mygenes_uniprot


# Code cell n°31

prot_universe <- unique(res_dHet_dHetRag_annot$gene_name)
prot_universe <- clusterProfiler::bitr(prot_universe,
                                       fromType = 'SYMBOL',
                                       toType = 'UNIPROT',
                                       OrgDb = 'org.Mm.eg.db')
prot_universe <- prot_universe[!duplicated(prot_universe$SYMBOL), ]
head(prot_universe)

# Code cell 32

enr_kegg <- enrichKEGG(gene = mygenes_uniprot$UNIPROT,
             keyType = "uniprot",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 1,
             universe = prot_universe$UNIPROT,
             organism = "mmu",
             pAdjustMethod = "none")

# Code cell n°33

dotplot(enr_kegg)

#Code cell n°34

options(repr.plot.width = 20, repr.plot.height = 15)
enr_kegg <- pairwise_termsim(enr_kegg, method = "JC", semData = "org.Mm.eg.db")
emapplot(enr_kegg, max.overlaps = 70, cex.params = list(category_label = 0.9))


# Code cell n°35

enr_kegg2 <- setReadable(enr_kegg, OrgDb = "org.Mm.eg.db", "UNIPROT")
cnetplot(enr_kegg2, categorySize="pvalue", showCategory = 5, cex.params = list(category_label = 1.1))


# Code cell n°36

sorted_genes <- res_dHet_dHetRag_annot[order(res_dHet_dHetRag_annot$log2FoldChange, decreasing = T),]

# To verify that the genes were sorted by decreasing Fold-change
head(sorted_genes)
tail(sorted_genes)

# We also remove duplicated genes 
sorted_genes <- sorted_genes[!duplicated(sorted_genes$gene_name), ]

# and we just save a vector of the orderedlog2FC with the gene_name as an attribute to the vector
sorted_genes_log2FC <- sorted_genes$log2FoldChange
names(sorted_genes_log2FC) <- sorted_genes$gene_name
str(sorted_genes_log2FC)
head(sorted_genes_log2FC)


# Code cell n°37

gsea_go <- gseGO(geneList = sorted_genes_log2FC, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.25, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none",
                 nPerm = 7500)    # by default there are 1000 permutations

# Code cell n°38

options(repr.plot.width = 15, repr.plot.height = 10) 
gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1])

# Code cell n°39

# this time we add the p-value for the enrichment of the pathway
gseaplot2(gsea_go, geneSetID = 2, title = gsea_go$Description[2], pvalue_table = TRUE)

# Code cell n°40

gseaplot2(gsea_go, geneSetID = 1:3, title = "top 3 gene sets", pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

# Code cell n°41

options(repr.plot.width = 17, repr.plot.height = 10)
dotplot(gsea_go, showCategory=20)+
        #label_format = function(x) stringr::str_wrap(x, width=  120))+
        ggtitle("dotplot for GSEA")

# Code cell n°42

options(repr.plot.width = 15, repr.plot.height = 10)

gsea_go <- enrichplot::pairwise_termsim(gsea_go, method = "JC", semData = "org.Mm.eg.db")
emapplot(gsea_go, showCategory = 20)

# Code cell n°43

options(repr.plot.width=10, repr.plot.height= 8)
    cnetplot(gsea_go, categorySize="pvalue", showCategory = 5)

# Code cell n°44

options(repr.plot.width = 18, repr.plot.height = 18)
cnetplot(gsea_go, categorySize="pvalue", showCategory = 5, max.overlaps = 100)

# Code cell n°45

options(repr.plot.width = 30, repr.plot.height = 10)
ridgeplot(gsea_go, label_format = function(x) stringr::str_wrap(x, width=  120))

# Code cell n°46

gsea_go2 <- gseGO(geneList = sorted_genes_log2FC, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.25, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none",
             nPermSimple = 7500, # by default there are 1000 permutations
             eps = 0)

# Code cell n°47

options(repr.plot.width = 15, repr.plot.height = 10) 
gseaplot2(gsea_go2, geneSetID = 1:3, title = "top 3 gene sets", pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

# Code cell n°48

options(repr.plot.width = 17, repr.plot.height = 15)
dotplot(gsea_go2, showCategory=20)+
        #label_format = function(x) stringr::str_wrap(x, width=  120))+
        ggtitle("dotplot for GSEA")

# Code cell n°49

options(repr.plot.width = 15, repr.plot.height = 10)

gsea_go2 <- enrichplot::pairwise_termsim(gsea_go2, method = "JC", semData = "org.Mm.eg.db")
emapplot(gsea_go2, showCategory = 20)

# Code cell n°50

options(repr.plot.width = 18, repr.plot.height = 18)
cnetplot(gsea_go2, categorySize="pvalue", showCategory = 5, max.overlaps = 100)

# Code cell n°51

options(repr.plot.width = 30, repr.plot.height = 15)
ridgeplot(gsea_go2, label_format = function(x) stringr::str_wrap(x, width=  120))

# Code cell n°52

str(norm_counts)

# Code cell n°53

gsea_norm_counts <- data.frame("NAMES" = norm_counts$gene_name, "DESCRIPTION" = NA, norm_counts[3:13])
head(gsea_norm_counts, n = 5)

dim(gsea_norm_counts)

str(gsea_norm_counts)
write.table(gsea_norm_counts, file=paste0(gseafolder, "gsea_norm_counts.txt"), sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# Code cell n°54

str(sorted_genes)
head(sorted_genes, n = 4)
tail(sorted_genes, n = 4)

# Code cell n°55

gsea_dHet_dHetRag <- data.frame(sorted_genes$gene_name, sorted_genes$log2FoldChange)
write.table(gsea_dHet_dHetRag, file=paste0(gseafolder, "gsea_dHet_dHetRag.rnk"), sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# Code cell n°56

str(samples)
samples

# Code cell n°57

gsea_pheno <- file(paste0(gseafolder, "gsea_pheno.cls"))
my_text <- paste0("11 3 1\n#dHet dHetRag WT\n", paste(samples$Condition, collapse=" "))
writeLines(my_text, gsea_pheno)
close(gsea_pheno)

# Code cell n°58

install.packages("flashClust")
install.packages("WGCNA")

library(flashClust)
library(WGCNA)

sessionInfo()

# Code cell n°59

enableWGCNAThreads(nThreads = 6) # to be adjusted to the exact size of your session

# Code cell n°60

hclust <- flashClust::hclust

# Code cell n°61

# This step is necessary to avoid an error arising from the use of different cor() functions.
# Fix found thanks to a japanese post :-) https://qiita.com/Razumall/items/8aa6417f78c4857b3670 
# More info there: https://programmersought.com/article/90752004413/

cor <- WGCNA::cor

# Code cell n°62
 
data <- norm_counts[1:13]

# You see that genes are listed in a column named "ensemblID" and samples are in columns
head(data) 

dim(data)

# Code cell n°63

datExpr = as.data.frame(t(data[, -c(1:2)]));
names(datExpr) = data$ensemblID;
rownames(datExpr) = names(data)[-c(1:2)];

# now samples are rows and genes are columns
head(datExpr, 13)

# Code cell n°64

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
 
 

# Code cell n°64b (optional)

#if (!gsg$allOK)
#   {if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
#       if (sum(!gsg$goodSamples)>0)
#           printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
#       datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
#       }
 

# Code cell n°65

sampleTree = hclust(dist(datExpr), method = "average");

options(repr.plot.width = 14, repr.plot.height = 7) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot <- plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
             cex = 1.5,
             cex.lab = 1.5, 
             cex.axis = 1.5, 
             cex.main = 2)

# Plot a line to show the cut: useless here as no outliers are present in our samples
abline(h = 160, col = "red")

# Code cell n°65b (optional)

# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

# Code cell n°66

#Create an object called "datTraits" that contains your trait data

datTraits0 <- read.table(paste0(gohome,"alldata/Data/sampleData-GSE158661-traits.tsv"),
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE)

head(datTraits0, 12)

# Code cell n°67

#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits0) <- datTraits0$SampleName
head(datTraits0, 12)

# Code cell n°68

datTraits0$SampleName = NULL
table(rownames(datTraits0)==rownames(datExpr))
head(datTraits0, 12)
 

# Code cell n°69

datTraits <- datTraits0[, -c(1:9)]
head(datTraits, 12)
dim(datTraits)

# Code cell n°70

save(datExpr, datTraits, file=paste0(wgcnafolder,"SamplesAndTraits.RData"))


# Code cell n°71

options(repr.plot.width = 14, repr.plot.height = 8) 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# Code cell n°72

collectGarbage()

# Code cell n°73 (optional)

# Load the expression and trait data saved in the first part
stdata <- paste0(wgcnafolder, "SamplesAndTraits.RData")
stdata
load(stdata, verbose = TRUE)

# Code cell n°74
 
# Choose a set of soft powers
powers = c(c(1:10), seq(from =10, to=37, by=1)) #choosing a set of soft-thresholding powers

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5) #call network topology analysis function

# For information, we could perform an analysis considering the network as "signed", but as we have few samples, the SFT index is reduced. 
# We would end with a higher value for the chosen power (36 instead of 33)
#sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

# Code cell n°75

options(repr.plot.width = 25, repr.plot.height = 10) 
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

collectGarbage()

# Code cell n°76

enableWGCNAThreads(nThreads = 6)
softPower <- 33
adjacency <- adjacency(datExpr, power = softPower, type = "unsigned") # specify network type
head(adjacency)
 
collectGarbage()

# Code cell n°77

#translate the adjacency into topological overlap matrix (TOM) and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="unsigned") # specify network type
dissTOM = 1-TOM
 
collectGarbage()
# Code cell n°78

save(dissTOM,softPower, file=paste0(wgcnafolder, "wgcna-dissTOM.RData"))

collectGarbage()# Optional Code cell 78bis **

file.copy(paste0(myhome, "alldata/Results/wgcna/wgcna-dissTOM.RData"), wgcnafolder))

# Code cell n°79

# if you want to start from the previously stored R objet, uncomment the following lines to load the dissTOM matrix
# and adjust the value used for power if needed
#ddata <- paste0(wgcnafolder, "wgcna-dissTOM.RData")
#ddata
#load(ddata, verbose = TRUE)
# softPower <- 33

# This tree plot is usually wide
options(repr.plot.width = 25, repr.plot.height = 15) 
 
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.05)

#This sets the minimum number of genes to cluster into a module
# here we have a large number of genes, so we set this parameter to a higher value (default = 30)
minModuleSize = 100

dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

# Code cell n°80

collectGarbage()
# Code cell n°76b (alternative)

net = blockwiseModules(datExpr, power = 30,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "wgcna-mouseTOM", 
                       verbose = 2)# Code cell n°77b (alternative)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
collectGarbage()
# Code cell n°81

save(dynamicMods, dynamicColors, geneTree, MEList, MEs, MEDiss, METree, file=paste0(wgcnafolder,"Network_allGenes_unsigned.RData"))
  

# Code cell n°82

options(repr.plot.width = 25, repr.plot.height = 10) 

plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

 

# Code cell n°83

#plot dendrogram with module colors below it

options(repr.plot.width = 30, repr.plot.height = 15) 


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels= FALSE, hang=0.04, 
                    addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

 

# Code cell n°84

save(MEs, moduleLabels, moduleColors, mergedColors, geneTree, file=paste0(wgcnafolder, "wgcna-mouse-modules.RData"))

# Code cell n°85

collectGarbage()
# Code cell n°86

# Load the expression and trait data saved in the first part
stdata <- paste0(wgcnafolder, "SamplesAndTraits.RData")
stdata
load(stdata, verbose = TRUE)


# Load network data saved in the second part.
wdata <- paste0(wgcnafolder, "wgcna-mouse-modules.RData")
wdata
load(wdata, verbose = TRUE)

ndata <- paste0(wgcnafolder, "Network_allGenes_unsigned.RData")
ndata
load(ndata, verbose = TRUE)

pdata <- paste0(wgcnafolder, "wgcna-plotTOM.RData")
pdata
load(pdata, verbose = TRUE)
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

nGenes
nSamples

# Code cell n°87

# if you want to start from the stored R objet, uncomment the following lines to load the dissTOM matrix
#ddata <- paste0(wgcnafolder, "wgcna-dissTOM.RData")
#ddata
#load(ddata, verbose = TRUE)

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

collectGarbage()
# Code cell n°88

save(plotTOM, file=paste0(wgcnafolder, "wgcna-plotTOM.RData"))
collectGarbage()# Optional Code cell 88bis **

file.copy(paste0(myhome, "alldata/Results/wgcna/wgcna-plotTOM.RData"), wgcnafolder))

# Code cell n°89

# if you want to start from the stored R objet, uncomment the following lines to load the plotTOM matrix
#pdata <- paste0(wgcnafolder, "wgcna-plotTOM.RData")
#pdata
#load(pdata, verbose = TRUE)

# Call the plot function
options(repr.plot.width = 15, repr.plot.height = 10) 
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# Code cell n°90

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]


# Taking the dissimilarity to a power, say 7 or 10, makes the plot more informative 
# by effectively changing the color palette; 
# Setting the diagonal to NA also improves the clarity of the plot
options(repr.plot.width = 10, repr.plot.height = 10)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# Code cell n°91

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$Weight_g)
names(weight) = "weight"

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))

# Plot the relationships between the eigengenes and the trait
options(repr.plot.width = 10, repr.plot.height = 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)

# Code cell n°92

# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", 
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)

# Code cell n°93

collectGarbage()
# Code cell n°94 (optional)

# Load the expression and trait data saved in the first part
stdata <- paste0(wgcnafolder, "SamplesAndTraits.RData")
stdata
load(stdata, verbose = TRUE)

# Load network data saved in the second part and matrix data saved in the third part
wdata <- paste0(wgcnafolder, "wgcna-mouse-modules.RData")
wdata
load(wdata, verbose = TRUE)

ndata <- paste0(wgcnafolder, "Network_allGenes_unsigned.RData")
ndata
load(ndata, verbose = TRUE)

ddata <- paste0(wgcnafolder, "wgcna-dissTOM.RData")
ddata
load(ddata, verbose = TRUE)

pdata <- paste0(wgcnafolder, "wgcna-plotTOM.RData")
pdata
load(pdata, verbose = TRUE)
# Code cell n°95

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

nGenes
nSamples

# Code cell n°96

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

collectGarbage()

# Code cell n°97

options(repr.plot.width = 25, repr.plot.height = 20) 
par(mfrow = c(1,2))

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))


# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(20),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Code cell n°98

# Define variables weight and survival containing the 2 columns of datTrait
weight = as.data.frame(datTraits$Weight_g);
names(weight) = "weight"

survival <- as.data.frame(datTraits$Survival_days);
names(survival) = "Survival"

# names (colors) of the modules computed above
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# Computing the significance of the correlation between the genes and a specific trait
geneTraitSignificanceW = as.data.frame(cor(datExpr, weight, use = "p"))
GSPvalueW = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceW), nSamples))

names(geneTraitSignificanceW) = paste("GS.W.", names(weight), sep="")
names(GSPvalueW) = paste("p.GS.W.", names(weight), sep="")

geneTraitSignificanceS = as.data.frame(cor(datExpr, survival, use = "p"))
GSPvalueS = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceS), nSamples))

names(geneTraitSignificanceS) = paste("GS.S.", names(survival), sep="")
names(GSPvalueS) = paste("p.GS.S.", names(survival), sep="")

# Code cell n°99

module = "lightgreen"
column = match(module, modNames)
moduleGenes = moduleColors==module


options(repr.plot.width = 20, repr.plot.height = 10) 

par(mfrow = c(1,2));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceW[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceS[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for survival",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Code cell n°100

names(datExpr)[moduleColors=="lightgreen"]

# Code cell n°101

table(norm_counts[1]==rownames(geneModuleMembership))

# Code cell n°102

#Create the starting data frame
geneInfo0 = data.frame(norm_counts[1],
                      geneSymbol = norm_counts[2],
                      moduleColor = moduleColors,
                      geneTraitSignificanceW,
                      GSPvalueW)

# Code cell n°103

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

head(geneInfo0)

# Code cell n°104

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.W.weight));
geneInfo = geneInfo0[geneOrder, ]

# Code cell n°105

write.table(geneInfo, file = paste0(wgcnafolder, "geneInfo-Weight.tsv"), sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# Code cell n°106

#Create the starting data frame
geneInfo0 = data.frame(norm_counts[1],
                      geneSymbol = norm_counts[2],
                      moduleColor = moduleColors,
                      geneTraitSignificanceS,
                      GSPvalueS)

# Order modules by their significance for survival
modOrder = order(-abs(cor(MEs, survival, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Code cell n°107

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.S.Survival));
geneInfo = geneInfo0[geneOrder, ]

# Code cell n°108

write.table(geneInfo, file = paste0(wgcnafolder, "geneInfo-Survival.tsv"), sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

## Code cell 109 ##   

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE)

runfolder <- paste0(myfolder,"/run_notebooks")
       
file.copy(paste0(myfolder, "/Pipe_11-R-ORA-GSEA.ipynb"), paste0(runfolder, "/Pipe_11-R-ORA-GSEA-run.ipynb"))

