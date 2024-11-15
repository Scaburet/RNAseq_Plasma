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

cat("Voici mon environnement de travail avec les paquets de R chargés:\n")
sessionInfo()

## Code cell 3 ##


gohome <- "/shared/projects/2413_rnaseq_cea/"
gohome

# In a Jupyter Hub and a jupyter notebook in R, by default the working directory is where the notebook is opened
getwd()
myfolder <- getwd()
myfolder


## Code cell 4 ##

# storing the path to this folder in a variable
rintrofolder <- paste0(myfolder,"/Results/Rintro/", sep = "")
rintrofolder

# listing the content of the folder
print(system(paste("ls -hlt", rintrofolder), intern = TRUE) )

## Code cell 5 ##

options(repr.plot.width=10, repr.plot.height=10)

## Code cell 6 ##

rdata <- paste0(rintrofolder,"myDataf.RData")
rdata
load(rdata,verbose = T)

## Code cell 7 ##

ls() 

## Code cell n°8 ##

plot(myDataf$weight ~ myDataf$size)  

## Code cell n°9 ##

boxplot(myDataf$weight)

## Code cell n°10 ##

boxplot(myDataf$weight ~ myDataf$sex) 

## Code cell n°11 ##

a <- rnorm(1000) # to sample 1000 values from a normal distribution of mean 0 and standard deviation 1
hist(a, breaks = 20) # the argument breaks is used to specify the number of intervals

## Code cell n°12 ##

plot(myDataf$weight ~ myDataf$size,   # this is a primary function
     main = "Weight ~Size", 
     xlim = c(-3,3),
     ylim = c(0,200),
     type = "n",
     xlab = "Size",
     ylab = "Weight")
points(myDataf$weight[1:2] ~ myDataf$size[1:2], # this is a secondary function
       pch = 6,                                 # this is a parameter of the points() function: pch = plot character, i.e shape of the points.  6 is a downward pointing triangle
       col = "blue")
points(myDataf[3:6,"weight"] ~ myDataf$size[3:6],
       type = "b",
       pch = 23,
       col = "magenta",
       bg = "cyan",
       cex = 2)
points(seq(0, 2.5, 0.5),
       c(1, 10, 25, 50, 125, 150),
       type = "l") 
lines(-seq(0, 2.5, 0.5),
      c(1, 10, 25, 50, 125, 150),
      lty = "dotdash",
      col = "blue",
      lwd = 3) 
abline(0, -50,
       lty = 3,
       col = "red")
abline(v = 0,
       lty = 2,
       col = "green")

mtext("overlap of unrelated graphs",
      side = 1)
mtext("other text",
      side = 1,
      line = 2)
text(-1, 150,
     "some text here")
axis(side = 4,
     labels = c(0, 20, 80, 160),
     at = c(0, 20, 80, 160),
     tick = TRUE)
legend("topright",
       c("blue triangles", "diamonds", "red line"),
       pch = c(6, 23, NA),
       col = c("blue", "magenta", "red"),
       pt.bg = c("transparent", "cyan", "transparent"),
       lty = c(0,0,3))

### Code cell n°13 ##

summary(myDataf)

## Code cell n°14 ##

t.test(myDataf$weight ~ myDataf$sex) 

## Code cell n°15 ##

cases_controls <- data.frame("genotype" = c("0/0", "1/0 (or 0/1)","1/1"),
                             "cases" = c(188,133,24),
                             "controls" = c(94,92,30))

## Code cell n°16 ##

cases_controls

## Code cell n°17 ##

proportions(as.matrix(cases_controls[,2:3]),2)

## Code cell n°18 ##

barplot(proportions(as.matrix(cases_controls[,2:3]),2),
        beside = FALSE,
        legend = cases_controls$genotype)

## Code cell n°19 ##

chisq.test(cases_controls[, 2:3]) 

## Code cell n°20 ##

motorisation <- read.table(paste0('/shared/projects/2413_rnaseq_cea/alldata/Example_Data/motorisation.txt'), header = FALSE, stringsAsFactors = FALSE) # nommez votre objet R
#raw cell
#Alternativement, si le fichier motorisation est dans votre repertoire de travail ou un autre, donnez son chemin absolu ou relatif
motorisation <- read.table("./motorisation.txt", header = F, stringsAsFactors = F) 
## Code cell n°21 ##
str(motorisation)

## Code cell n°22 ##
names(motorisation) <- "type_de_motorisation"

## Code cell n°23 ##
table(motorisation$type_de_motorisation)

## Code cell n°23.2 ##
# Si besoin, vous pouvez extraire les noms attribués à ces valeurs avec la commande
names(table(motorisation$type_de_motorisation))

## Code cell n°24 ##
pie(table(motorisation$type_de_motorisation))

## Code cell n°25 ##
pie(table(motorisation$type_de_motorisation), col = c("green3","blue","magenta","orange"))

## Code cell n°25.2 ##
?pie

## Code cell n°25.3 ##
pie(c(9, 1, 7, 5), labels = c("Diesel","Electrique", "Essence", "Hybride"),
    col = c("green3","blue","magenta","orange"))

## Code cell n°25.4 ##
pie(sort(table(motorisation$type_de_motorisation)), col = c("green3","blue","magenta","orange")) # les couleurs sont attribuées dans le nouvel ordre...qui est peut etre mieux car diesel n'est plus en vert!

## Code cell n°26 ##
barplot(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation),
        col = c("green3","blue","magenta","orange"))

## Code cell n°26.2 ##
barplot(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation),
        cex.names = 1.5,
        col = c("green3","blue","magenta","orange"))

## Code cell n°26.3 ##
barplot(matrix(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation)),
        col = c("green3","blue","magenta","orange"),
        legend.text = unique(sort(motorisation$type_de_motorisation)),
        args.legend = list(x = "topright", cex = 0.8, bty = "n"),
        width = 1,
        xlim = c(0,5) )

## Code cell n°26.4 ##
barplot(matrix(sort(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation),decreasing = T)),
        col = c("green3","blue","magenta","orange"),
        legend.text = unique(sort(motorisation$type_de_motorisation)),
        args.legend = list(x = "topright", cex = 0.8, bty = "n"),
        width = 1,
        xlim = c(0,5) )

## Code cell n°27 ##
opar <- par(no.readonly=TRUE) # l'argument no.readonly=TRUE permet de supprimer l'affichage d'éventuels warnings 
par(mfrow = c(1,2)) 
pie(table(motorisation$type_de_motorisation), col = c("green3","blue","magenta","orange"))
barplot(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation),col = c("green3","blue","magenta","orange"))


## Code cell n°27.2 ##
opar <- par(no.readonly=TRUE)
par(mfrow = c(1,2)) 
pie(sort(table(motorisation$type_de_motorisation)),
    col = c("green3","blue","magenta","orange"),
    clockwise = T)  
barplot(matrix(sort(table(motorisation$type_de_motorisation)/length(motorisation$type_de_motorisation),
                     decreasing = T)),
         col = c("orange","magenta","blue","green3"),
         legend.text = names(sort(table(motorisation$type_de_motorisation)
                                  ,decreasing = T)),
         args.legend = list(x = "topright", cex = 0.8, bty = "n"),
         width = 1,
         xlim = c(0,5))


## Code cell n°28 ##
myrandomdata <- rnorm(100, mean = 10, sd = 5)# bien penser à assigner le résultat de votre tirage dans un objet, sinon, les valeurs changent à chaque fois que vous effectuez un tirage avec la commande rnorm!

## Code cell n°29 ##

## Code cell n°30 ##

## Code cell n°31 ##

## Code cell n°32 ##
opar <- par()
par(mfrow = c(4,1))
hist(myrandomdata, breaks = 5)
hist(myrandomdata, breaks = 50)
hist(myrandomdata, breaks = 100)
boxplot(myrandomdata, horizontal = T)
suppressWarnings(par(opar))

## Code cell n°33 ##
round(pnorm(q = 15, mean = 10, sd = 5) - pnorm(q = 7, mean = 10, sd = 5), 3)  

## Code cell n°34 ##
round(1 - pnorm(q = 15, mean = 10, sd = 5, lower.tail = FALSE) - pnorm(q = 7, mean = 10, sd = 5), 3)

## Code cell n°35 ##
round(qnorm (p = 0.67, mean = 10, sd = 5), 1)

## Code cell n°36 ##
round(pchisq(q = 6.26, 3, lower.tail = F), 4)

## Code cell n°37 ##
round(pchisq(3.84,1, lower.tail=F), 4)

## Code cell n°38 ##
qchisq(0.05, 1, lower.tail=F)

## Code cell n°39 ##
data(airquality)
str(airquality)

## Code cell n°40 ##
opar <- par()
par(mfrow = c(2,2),
    mar = c(4.1, 2.1, 5.1, 2.1))
stripchart(airquality$Temp~airquality$Month,
           pch = 22,
           xlab = "Temperature", ylab = "Month",
           main = "A. stripchart of temperatures",
           col = 6)
hist(airquality$Temp,
     xlab = "Temp",
     freq = FALSE,
     col="transparent",
     main = "B. histogram of temperatures")
lines(density(airquality$Tem), col = 5)
boxplot(airquality$Ozone~airquality$Month,
        names = c("May", "June", "July", "August", "September"),
        cex.axis = 0.7,
        col = 4,
        pch = "'",
        main = "C. boxplot of ozon level per month")
plot(airquality$Ozone~airquality$Wind,
     xlab = "Wind",
     ylab = "Ozone",
     xlim = c(5,14),
     pch = "'",
     main = "D. ozone level versus wind")
abline(lm(airquality$Ozone~airquality$Wind),
       col = 2, lty = 2, lwd = 1.5)
title("Graphs from the airquality dataset", outer = T, line = -1)
suppressWarnings(par(opar))

## Code cell n°40.2 ##
fmonth <- factor(airquality$Month,levels = 5:9)
levels(fmonth) <- c("May","June","July","August","September")
boxplot(airquality$Ozone ~ fmonth ,
        col = "blue",
        ylab = "Ozone",
        main = "C. boxplot of ozon level per month")

## Code cell n°40.3 ##
boxplot(data = airquality, Ozone~fmonth , col = "blue", ylab = "Ozone", main = "C. boxplot of ozon level per month")

## Code cell n°40.4 ##
attach(airquality)
boxplot(Ozone ~ fmonth , col = "blue", ylab = "Ozone", main = "C. boxplot of ozon level per month")
detach(airquality)

## Code cell n°41 ##
sample1 <- rnorm(n = 100, mean = 174, sd = 7)

## Code cell n°42 ##
sample2 <- rnorm(n = 100, mean = 162, sd = 7)

## Code cell n°43 ##
summary(sample1)
summary(sample2)

## Code cell n°44 ##
boxplot(sample1, sample2)

## Code cell n°45 ##
t.test(sample1, sample2, alternative = "two.sided")

## Code cell n°46 ##
str(t.test(sample1, sample2, alternative = "two.sided"))

## Code cell n°46.2 ##
t.test(sample1, sample2, alternative = "two.sided")$p.value

## Code cell n°47 ##
compute_tv <- function(n1, m1, s1, n2, m2, s2){
    sample1 <- rnorm(n = n1, mean = m1, sd = s1)
    sample2 <- rnorm(n = n2, mean = m2, sd = s2)
    tval <- t.test(sample1, sample2, alternative = "two.sided")$statistic
    pval <- t.test(sample1, sample2, alternative = "two.sided")$p.value
    tv <- list("t" = tval, "p" = pval)
    return(tv)
}

## Code cell n°48 ##
compute_tv(n1 = 100, m1 = 174, s1 = 7.1, n2 = 100, m2 = 162, s2 = 6.5)

## Code cell n°49 ##
compute_tv(n1 = 30, m1 = 174, s1 = 7.1, n2 = 30, m2 = 162, s2 = 6.5)

## Code cell n°50 ##
compute_tv(30, 174, 25, 30, 162, 25)

## Code cell n°51 ##
compute_tv(30, 163, 3, 30, 160, 3)

## Code cell n°52 ##
df <- data.frame(height = c(sample1, sample2), grp = rep(c(1,2), each = 100))

## Code cell n°53 ##
head(df)

## Code cell n°54 ##
tail(df)

## Code cell n°55 ##
str(df)

## Code cell n°56 ##
table(df$grp)

## Code cell n°57 ##
tvalues <- c()
pvalues <- c()

for(i in 1:10000){
    perm_df <- df
    perm_df$grp <- sample(perm_df$grp, size=200)
    tval <- t.test(perm_df$height ~ perm_df$grp, alternative = "two.sided")$statistic
    tvalues <- c(tvalues, tval)
    pval <- t.test(perm_df$height ~ perm_df$grp, alternative = "two.sided")$p.value
    pvalues <- c(pvalues, pval)
}

## Code cell n°58 ##
summary(tvalues)

## Code cell n°59 ##
hist(tvalues)

## Code cell n°60 ##
qnorm(0.025, 0, 1, lower.tail=FALSE)

## Code cell n°61 ##
length(which(abs(tvalues)>=1.96))

## Code cell n°62 ##
summary(pvalues)
hist(pvalues)

## Code cell n°63 ##
coregone <- read.table("/shared/projects/2413_rnaseq_cea/alldata/Example_Data/poisson.txt", stringsAsFactors = FALSE, header = TRUE)
str(coregone)

## Code cell n°64 ##
summary(coregone)

## Code cell n°65 ##
cat("moyenne des variables:\n")
apply(coregone, 2, mean)# apply est une fonction qui permet d'appliquer des fonctions par colonne ou ligne d'un dataframe

## Code cell n°66 ##
cat("Ecart type des variables:\n")
apply(coregone, 2, sd)

## Code cell n°67 ##
##histogrammes avec densite des donnees en rouge et loi normale supperposee en cyan
opar <- par()
par(mfrow = c(3, 2))

hist(coregone$longueur_mm,prob=T,
     main = "Histogramme de la longueur totale",
     xlab = colnames(coregone)[1],
     ylab = "Densit?")
lines(density(coregone$longueur_mm),col="red" )
curve(dnorm(x,mean(coregone$longueur_mm), sd(coregone$longueur_mm)), col="cyan", add=T)

hist(coregone$poids_g,prob=T,
     main = "Histogramme du poids",
     xlab = colnames(coregone)[2],
     ylab = "Densit?", ylim=c(0,0.02))
lines(density(coregone$poids_g), col="red" )
curve(dnorm(x,mean(coregone$poids_g),sd(coregone$poids_g)), col="cyan", add=T)

hist(coregone$poids_gonades_mg,prob=T,
     main = "Histogramme du nombre de gonades",
     xlab = colnames(coregone)[3],
     ylab = "Densit?", ylim=c(0,0.02))
lines(density(coregone$poids_gonades_mg),col="red" )
curve(dnorm(x,mean(coregone$poids_gonades_mg), sd(coregone$poids_gonades_mg)), col="cyan",add=T)

hist(coregone$age_annee,prob=T,
     main = "Histogramme de l'age",
     xlab = colnames(coregone)[4],
     ylab = "Densit?", ylim=c(0,0.4))
lines(density(coregone$age_annee),col="red" )
curve(dnorm(x,mean(coregone$age_annee), sd(coregone$age_annee)), col="cyan",add=T)

hist(coregone$oeufs_nombre,prob=T,
     main = "Histogramme du nombre d'oeufs",
     xlab = colnames(coregone)[5],
     ylab = "Densit?", ylim=c(0,1.4e-04))
lines(density(coregone$oeufs_nombre),col="red" )
curve(dnorm(x,mean(coregone$oeufs_nombre), sd(coregone$oeufs_nombre)),col="cyan",add=T)
suppressWarnings(par(opar))

## Code cell n°68 ##
dnorm(mean(coregone$oeufs_nombre),mean(coregone$oeufs_nombre),sd(coregone$oeufs_nombre))

## Code cell n°69 ##
str(hist(coregone$oeufs_nombre, plot=F))

## Code cell n°69.2 ##
max(hist(coregone$oeufs_nombre, plot=F)$density) # la valeur la plus grande sur l'axe des Y

## Code cell n°70 ##
opar <- par()
par(mfrow = c(3, 2))
qqnorm(coregone$longueur_mm,main="QQ Plot Longueur")
qqline(coregone$longueur_mm)
qqnorm(coregone$poids_g,main="QQ Plot Poids")
qqline(coregone$poids_g)
qqnorm(coregone$poids_gonades_mg,main="QQ Plot Gonades")
qqline(coregone$poids_gonades_mg)
qqnorm(coregone$age_annee,main="QQ Plot Age")
qqline(coregone$age_annee)
qqnorm(coregone$oeufs_nombre,main="QQ Plot Oeufs")
qqline(coregone$oeufs_nombre)            
suppressWarnings(par(opar))

## Code cell n°71 ##
shapiro.test(coregone$longueur_mm)

## Code cell n°72 ##
shapiro.test(coregone$poids_g)

## Code cell n°73 ##
shapiro.test(coregone$poids_gonades_mg)

## Code cell n°74 ##
shapiro.test(coregone$age_annee)

## Code cell n°75 ##
shapiro.test(coregone$oeufs_nombre)

## Code cell n°76 ##
coregone$classe_age <- NA
coregone[which(coregone[,4] <= 10),"classe_age"] <- "jeune"
coregone[which(coregone[,4] > 10),"classe_age"] <- "vieux" 

## Code cell n°76.2 ##
table(coregone$classe_age)

## Code cell n°77##
boxplot(coregone$oeufs_nombre ~ coregone$classe_age) 
stripchart(coregone$oeufs_nombre ~ coregone$classe_age, col = "blue", add = TRUE, vertical = TRUE, pch = 20) 

## Code cell n°77.2 ##
print(boxplot(coregone$oeufs_nombre~coregone$classe_age, plot=F)$stats)

## Code cell n°77.3 ##
opar <- par()
par(mfrow=c(2,2))
boxplot(coregone$longueur_mm~coregone$age_annee) 
boxplot(coregone$poids_g~coregone$age_annee)
boxplot(coregone$poids_gonades_mg~coregone$age_annee)
boxplot(coregone$oeufs_nombre~coregone$age_annee)
suppressWarnings(par(opar))

## Code cell n°78 ##
tapply(coregone$oeufs_nombre, coregone$classe_age, mean, na.rm = TRUE)

## Code cell n°78.2 ##
# ou par vecteur en deux lignes:
mean(coregone$oeufs_nombre[which(coregone$classe_age == "jeune")])
mean(subset(coregone, classe_age == "vieux")$oeufs_nombre) # ou ici avec subset au lieu de which

## Code cell n°79 ##
t.test(coregone$oeufs_nombre~coregone$classe_age)

## Code cell n°80 ##
wilcox.test(coregone$oeufs_nombre~coregone$classe_age)

## Code cell n°81 ##
plot(coregone$oeufs_nombre ~ coregone$poids_gonades_mg, xlab = "Poids des gonades", ylab = "Nombre d'oeufs",pch = 20, main = "Relation entre le poids des gonades et le nombre d'oeufs")
abline(lm(coregone$oeufs_nombre ~ coregone$poids_gonades_mg),col="red")

## Code cell n°82 ##
cor.test(coregone$oeufs_nombre, coregone$poids_gonades_mg)

## Code cell n°83 ##
summary(lm(coregone$oeufs_nombre ~ coregone$poids_gonades_mg))

## Code cell 84 ##   

myfolder <- getwd()

#myfolder <- setwd('/shared/ifbstor1/projects/2413_rnaseq_cea/mylogin') # devrait être inutile
#myfolder

file.copy("/shared/projects/2413_rnaseq_cea/staff/Pipe_08-R_counts-exploratory-analysis-I.ipynb", myfolder)
## Code cell n°85 ##

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
#dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE)

runfolder <- paste0(myfolder,"/run_notebooks")
runfolder

file.copy(paste0(myfolder, "/Pipe_07b-R_plots_stats_withR.ipynb"), paste0(runfolder, "/Pipe_07b-R_plots_stats_withR_run.ipynb"))
