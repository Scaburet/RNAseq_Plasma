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
    "ggfortify",
    "ggrepel",
    "RColorBrewer",
    "ggplot2",
    "stringr",
    "matrixStats",
    "corrplot",
    "BiocManager",
    "FactoMineR",
    "factoextra"
)

# list the required libraries from the Bioconductor project
requiredBiocLib <- c(
    "DESeq2",
    "ComplexHeatmap")


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
# we can skip this step as the folder is already created
# dir.create(paste(myfolder,"/Results/deseq2/", sep = ""), recursive = TRUE)

# storing the path to this output folder in a variable
deseq2folder <- paste(myfolder,"/Results/deseq2/", sep = "")
deseq2folder

# listing the content of the folder
print(system(paste("ls -hlt", deseq2folder), intern = TRUE) )

## Code cell 7 ##

options(repr.plot.width = 15, repr.plot.height = 8) # for figure display in the notebook

## Code cell 8 ##

rdata <- paste0(deseq2folder, "deseq2.RData") # generate the path of the deseq2 data
rdata # check the path
load(rdata, verbose = TRUE) # load the path and ask to write what is loaded
rm(rdata) # remove the path of deseq2 data, no longer needed

## Code cell 9 ##

ls() 

samples

## Code cell 10 ##

head(rlog.dds2.annot) 
str(rlog.dds2.annot)

## Code cell 11 ##

head(rlog.dds2.annot[ , 19])
head(rlog.dds2.annot[,"gene_name"]) #the two commands are equivalent since column 19 name is "gene_name"

## Code cell 12 ##

length(unique(rlog.dds2.annot$ensemblID))
length(unique(rlog.dds2.annot$gene_name))
table(table(rlog.dds2.annot$gene_name)) # contingency table of gene_wise contingency tables
dup_genes <- names(which(table(rlog.dds2.annot$gene_name) == 2)) # to get the gene_names of the 16 gene_names with two ensemblID
subset(rlog.dds2.annot, gene_name %in% dup_genes)[, c("ensemblID", "gene_name")]

## Code cell 13 ##

#norm_counts <- rlog.dds2.annot[,c(19, 2:12)]
# or sidem in a more explicit way:
norm_counts <- rlog.dds2.annot[,c("ensemblID", "gene_name",
                                  "dHet_B-ALL_686_rep1", "dHet_B-ALL_686_rep2",
                                  "dHet_B-ALL_713_rep1", "dHet_B-ALL_713_rep2",
                                 "dHet_B-ALL_760_rep1", "dHet_B-ALL_760_rep2",
                                 "dHet_FetalLiver_proB_rep1",  "dHet_FetalLiver_proB_rep2", "dHet_FetalLiver_proB_rep3",
                                  "wt_BoneMar_proB_rep1", "wt_BoneMar_proB_rep2")]
dim(norm_counts)
head(norm_counts, n = 5)
summary(norm_counts[, -c(1:2)]) # useless to do summary on the first two columns that contains qualitative data

## Code cell 14 ##

rm(rlog.dds2.annot)
ls()

## Code cell 15 ##

# make a colour vector
conditionColor <- match(samples$Condition, c("dHet", "dHetRag", "WT")) + 1
# '+1' to avoid color '1' i.e. black

# Check distributions of samples using boxplots, using only the columns with read counts
boxplot(norm_counts[, 3:13], # or can do -c(1:2)
        xlab = "",
        ylab = "rlog.dds2.annot Counts",
        las = 2,
        col = conditionColor,
        main = "rlog.dds2.annot Counts")
# Let's add a blue horizontal line that corresponds to the median
abline(h = median(as.matrix(norm_counts[ , -c(1:2)])), col = "blue")

## Code cell 16 ##

# run PCA
PCAdata <- prcomp(t(norm_counts[, -c(1:2)])) # we get rid of the first tow columns of norm_counts that do not contain norm data
summary(PCAdata)

## Code cell 17 ##

# to display the two scree plots side by side
layout(matrix(1:2, ncol = 2))

screeplot(PCAdata) # barplot representation
screeplot(PCAdata, type = "lines", main = "Screeplot PCAdata - Eigenvalues") # same but with a line

## Code cell 18 ##

autoplot(PCAdata,
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
        geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), box.padding = 0.8)


## Code cell 19 ##

autoplot(PCAdata,
         x = 3,    # PC3
         y = 4,    # PC4
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
    geom_text_repel(aes(x = PC3, y = PC4, label = SampleName), box.padding = 0.8)


## Code cell 20 ##

length(unique(norm_counts$ensemblID))
## (Code cell 21) ##
# alternative way of finding the most variant genes with the coefficient of variation

# we define a function to compute the coefficient of variation
compute_cv <- function(x){
    sd(x,  na.rm=TRUE) / mean(x, na.rm=TRUE)
    }
    
# we apply our function to our matrix, and we store the results in cv_genes
var_genes <- apply(norm_counts[, 3:13], MARGIN = 1, FUN = compute_cv)
## Code cell 22 ##

var_genes <- matrixStats::rowVars(as.matrix(norm_counts[ , 3:13]))
#var_genes <- apply(X = norm_counts[, 3:13], MARGIN = 1, FUN = var, na.rm = TRUE) # simmilar way to generate var_genes

head(var_genes)
str(var_genes)
summary(var_genes)

## Code cell 23 ##

norm_counts$var <- var_genes
head(norm_counts)

## Code cell 24 ##

norm_counts <- norm_counts[order(norm_counts$var, decreasing = TRUE),]
norm_counts$rank_var <- 1:nrow(norm_counts)
head(norm_counts)

## Code cell 25 ##

top50var <- norm_counts[1:50,]
head(top50var)

## Code cell 26 ##

str(top50var)
dim(top50var)
length(unique(top50var$ensemblID))

head(top50var, n = 5)

## Code cell 27 ##

length(unique(top50var$gene_name))
top50var <- top50var[order(top50var$gene_name, decreasing = TRUE),]
head(top50var)

## Code cell 28 ##

top50var$gene_name

## Code cell 29 ##

row.names(top50var) <- top50var$gene_name

## Code cell 30 ##

str(res2_dHet_dHetRag_sig_ranked_annot)

## Code cell 31 ##

res2_dHet_dHetRag_sig_ranked_annot$rank_DE <- 1:nrow(res2_dHet_dHetRag_sig_ranked_annot)
str(res2_dHet_dHetRag_sig_ranked_annot)

## Code cell 32 ##

merge(top50var, res2_dHet_dHetRag_sig_ranked_annot[,-14],
      by = "ensemblID", all = FALSE, sort = FALSE) [,c("ensemblID", "gene_name", "rank_var", "rank_DE", "log2FoldChange", "padj")]

## Code cell 33 ##

setdiff(top50var$gene_name, res2_dHet_dHetRag_sig_ranked_annot$gene_name)

## Code cell 34 ##

top50DE <- res2_dHet_dHetRag_sig_ranked_annot[1:50, ]
setdiff(top50DE$gene_name, top50var$gene_name) # list the top50 DE genes that are not among the top50 var genes
merge(top50DE, norm_counts[, - 2], by = "ensemblID", all.x = TRUE, all.Y = TRUE, sort = FALSE)[,c("ensemblID", "gene_name", "rank_var", "rank_DE", "log2FoldChange", "padj")]

## Code cell 35 ##

# run PCA
PCAdata2 <- prcomp(t(top50var[, 3:13]))


## Code cell 36 ##

# to display the two scree plots side by side
layout(matrix(1:2, ncol = 2))

screeplot(PCAdata2)
screeplot(PCAdata2, type = "lines")

## Code cell 37 ##

autoplot(PCAdata2,
         data = samples, 
         colour = "Condition", 
         shape = "Tissue",
         size = 6) +
        geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), box.padding = 0.8)


## Code cell 38 ##

biplot(PCAdata2,
       scale = 0)

## Code cell 39 ##

for_factominer <- t(top50var [, 3:13])
str(for_factominer)
head(for_factominer)

## Code cell 40 ##

row.names(for_factominer) == samples$SampleName
for_factominer <- data.frame(for_factominer, samples)
str(for_factominer)

## Code cell 41 ##

for_factominer$mouseID <- c(rep("686", 2), rep("713", 2), rep("760", 2), 1:5)
for_factominer$B_type <- c(rep("B-ALL", 6), rep("proB", 5))

## Code cell 42 ##

summary(for_factominer$Xist)
hist(for_factominer$Xist)
for_factominer$sex <- ifelse(for_factominer$Xist < 6, "male", "female")
table(for_factominer$sex)
for_factominer[, c("Xist", "sex")]

## Code cell 43 ##

for_factominer$Sos1 <- unlist(subset(norm_counts, gene_name == "Sos1")[,3:13])
for_factominer$Gm2629 <- unlist(subset(norm_counts, gene_name == "Gm2629")[,3:13])
head(for_factominer)
str(for_factominer)

## Code cell 44 ##

# Use PCA() of FactoMineR
#?PCA
PCAres <- FactoMineR::PCA(for_factominer,
                          quali.sup = 51:(ncol(for_factominer)-2), # specify the index of the qulatitative colums
                          quanti.sup = (ncol(for_factominer)-1):ncol(for_factominer),
                                                                graph = FALSE)
#str(PCAres) # longue structure!
print(class(PCAres))
names(PCAres)

## Code cell 45 ##

# Creates graphs with the samples (called "individuals in factominer) according to each axsis :

FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", cex = 1)

## Code cell 46 ##

# Same graph adding colors for qualitative variables
# save in pdf if wanted

#pdf("PCA_individus.pdf")

FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", cex = 1)
FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", habillage = "Condition", cex = 1)
FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", habillage = "Tissue", cex = 1)
FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", habillage = "mouseID" , cex = 1)
FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", habillage = "B_type" , cex = 1)
FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", habillage = "sex", cex = 1)
#dev.off()

## Code cell 47 ##

FactoMineR::plot.PCA(PCAres, axes = c(2,3), choix = "ind", autoLab = "yes", invisible = "quali", cex = 1)

## Code cell 48 ##

PCAres$eig

# only for first five:
PCAres$eig[1:5]

## Code cell 49 ##

# Graphic Inertia and dimensions

#pdf("PCA_inertia.pdf")
eig.val <- PCAres$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        ylim = c(0,60),
        col ="steelblue")

barplot(eig.val[, 3], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Cumulative Percentage of variance",
        ylim = c(0,100),
        col ="steelblue")

#dev.off()

## Code cell 50 ##

# diagramme des éboulis avec 10 composantes par défaut
factoextra::fviz_eig(PCAres, addlabels = TRUE)

## Code cell 51 ##

# idem same with only 5 components:
factoextra::fviz_eig(PCAres, addlabels = TRUE, ncp = 5)

## Code cell 52 ##

round(PCAres$quali.sup$eta2, 2)

## Code cell 53 ##

round(PCAres$quanti.sup$cor,2)

## Code cell 54 ##

#pdf("nom_fichier.pdf")
# les composantes par défaut sont les 1 et 2
# vous pouvez les modifier via l'argument "axes"

FactoMineR::plot.PCA(PCAres, choix = "ind", autoLab = "yes", invisible = "quali", cex = 1,
                     habillage = "Condition", axes = c(1,2))

#dev.off()

## Code cell 55 ##

factoextra::fviz_pca_ind(PCAres, label = "none",
                         habillage = as.factor(PCAres$call$quali.sup$quali.sup$Condition),
             addEllipses = TRUE, ellipse.level = 0.95)

## Code cell 56 ##

factoextra::fviz_pca_ind(PCAres, label = "none",
                         habillage = as.factor(PCAres$call$quali.sup$quali.sup$sex),
             addEllipses = TRUE, ellipse.level = 0.95)

factoextra::fviz_pca_ind(PCAres, label = "none",
                         habillage = as.factor(PCAres$call$quali.sup$quali.sup$sex),
             addEllipses = TRUE, ellipse.level = 0.95, axes = c(2,3))

## Code cell 57 ##

factoextra::fviz_pca_ind(PCAres, label = "none",
                         habillage = as.factor(PCAres$call$quali.sup$quali.sup$Tissue),
             addEllipses = TRUE, ellipse.level = 0.95)

## Code cell 58 ##

factoextra::fviz_pca_ind(PCAres, label="none",
                         col.ind = as.numeric(PCAres$call$quanti.sup$Sos1))

## Code cell 59 ##

factoextra::fviz_pca_ind(PCAres, label="none",
                         col.ind = as.numeric(PCAres$call$quanti.sup$Gm2629))

## Code cell 60 ##

factoextra::fviz_contrib(PCAres, choice = "ind", axes = 1:2, top = 30)

## Code cell 61 ##

# Contributions to axis one of the top 10 genes contributing to it :
factoextra::fviz_contrib(PCAres, choice = "var", axes = 1, top = 10)

# Contributions to second axis:
factoextra::fviz_contrib(PCAres, choice = "var", axes = 2, top = 10)

## Code cell 62 ##

# Correlation circle
# with the argument "contrib 10" you display the top 10 genes with the best correlations with axes 1 and 2

#pdf(file = "PCA_diab_contrib30.pdf")
FactoMineR::plot.PCA(PCAres, choix = "var", cex = 1, select = "contrib 10", unselect = 1)
#dev.off()

## Code cell 63 ##

factoextra::fviz_pca_var(PCAres, select.var = list(contrib = 10))

## Code cell 64 ##

factoextra::fviz_pca_var(PCAres, select.var = list(cos2 = 0.8))

## Code cell 65 ##

factoextra::fviz_pca_biplot(PCAres, select.var = list(contrib = 10),
                            col.var = "black",
                            habillage = as.factor(PCAres$call$quali.sup$quali.sup$Condition))

## Code cell 66 ##

# Lignes de codes pour calculer les coefficients de corrélation entre le niveau des transcrits et la position des échantillons sur les axes
cor.test(for_factominer$Xist, PCAres$ind$coord[,1])
cor.test(for_factominer$Xist, PCAres$ind$coord[,2])

## Code cell 67 ##
   
clusters <- hclust(dist(as.matrix(t(top50var[ ,3:13]))), method = "ward.D" )
plot(clusters)

rm(clusters)    


## Code cell 68 ##

heatmap(as.matrix(top50var[,c(3:13)])) 

## Code cell 69 ##

t_top50var <- t(apply(top50var[,c(3:13)], 1, scale))

## Code cell 70 ##

# we change the dimension of the output for better rendering
options(repr.plot.width = 10, repr.plot.height = 10)

ComplexHeatmap::Heatmap(t_top50var,
                        name = "Z-score",
                        clustering_distance_rows = "pearson",
                        column_title = "pre-defined distance method (1 - pearson)", row_title = "Top 50 genes",
                        clustering_method_rows = "ward.D")

## Code cell 71 ##

heatmap(cor(t(top50var[,c(3:13)])))

## Code cell 72 ##

options(repr.plot.width = 13, repr.plot.height = 13)

# A first "simple" plot
corrplot(cor(t(top50var[,c(3:13)])))

# method = "square" is the default 
# hclust clusterises the genes, grouping together the most correlated ones
# rectangles around clusters and correlation coefficients are also added
corrplot(cor(t(top50var[,c(3:13)])),
         method = "square",
         order = 'hclust',
         addrect = 4,
         addCoef.col = 'black',
         number.cex = 0.3)

## Code cell 73 ##

# type = 'upper' will show only the upper half of the matrix
corrplot(cor(t(top50var[,c(3:13)])),
        method = 'ellipse',
        type = 'upper',
        insig = 'blank')

## Code cell 74 ##

print(ls())

## Code cell 74 ##

head(res2_dHet_dHetRag_sig_ranked_annot)

## Code cell 75 ##

rm(PCAdata, PCAdata2, t_top50var)

## Code cell 76 ##

ls()
save.image(file=paste0(deseq2folder,"deseq2-final.RData"))

## Code cell 77 ##   

myfolder

file.copy("/shared/projects/2413_rnaseq_cea/pipeline/Pipe_11-R-ORA-GSEA.ipynb", myfolder)


## Code cell 78 ##   

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE, showWarnings = FALSE)

runfolder <- paste0(myfolder,"/run_notebooks")
       
file.copy(paste0(myfolder, "/Pipe_10-R-Normcounts-exploratory-analysis-II.ipynb"), paste0(runfolder, "/Pipe_10-R4-Normcounts-exploratory-analysis-II-run-run.ipynb"))

