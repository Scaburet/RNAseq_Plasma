## Code cell n°1 ##

2+2

## Code cell n°2 ##

2-3

## Code cell n°3 ##

6/3

## Code cell n°4 ##

10/3

## Code cell n°5 ##

10%%3

## Code cell n°6 ##

mean(c(1,2)) # we will see further down in this notebook that we need to put concatenate different values with a c() first
exp(-2)

## Code cell n°7 ##

round(exp(-2), 2)

## Code cell n°8 ##

log(100, base = 10) # we want to get the log of 100 in base 10 

## Code cell n°9 ##

help(round)

## Code cell n°10 ##

?exp

## Code cell n°11 ##

x <- 2

## Code cell n°12 ##

x

## Code cell n°13 ##

x + x

## Code cell n°14 ##

y <- x + 3

## Code cell n°15 ##

y

## Code cell n°16 ##

y <- x + 3
y

## Code cell n°17 ##

s <- "this is a string of characters"
s

## Code cell n°18 ##

class(x)
class(s)

## Code cell n°19 ##

"1"
class("1")
class(1)

## Code cell n°20 ##

try("1" + 3)# I added the try function to avoid stopping the notebook if you want to run all the cells

## Code cell n°21 ##

1 + 3

## Code cell n°22 ##

sessionInfo()

## Code cell n°23 ##

getwd()

## Code cell n°24 ##

initial_path <- getwd()
initial_path #current working directory

## Code cell n°25 ##

setwd(paste0(initial_path,"/Results/", sep = "")) #change current working directory
getwd() #change is visible

## Code cell n°26 ##

setwd(initial_path) # reset to initial working directory by using the content of the variable
rm(initial_path)
getwd()

## Code cell n°27 ##

ls()

## Code cell n°28 ##

rm(y) # the object y is removed
ls()

## Code cell n°29 ##

dir()

## Code cell n°30 ##

list.files(pattern=".ipynb")

## Code cell n°31 ##

save(x,file = "x.RData")

## Code cell n°32 ##

rm(x)
ls()

## Code cell n°33 ##

load("x.RData")
ls()
x    # x is again accessible

## Code cell n°34 ##

file.remove("x.RData") #remove file: returns TRUE on successful removal

## Code cell n°35 ##

save(x, s, file = "xands.RData")

## Code cell n°36 ##

file.remove("xands.RData")  # to clean the working directory

## Code cell n°37 ##

ls()
save.image(file = "AllMyData.RData")

## Code cell n°38 ##

rm(list=ls()) # this command removes all the objects on the R session
ls() #all variables have been removed

## Code cell n°39 ##

load("AllMyData.RData")
ls() #all variables are accessible again
file.remove("AllMyData.RData")
ls()
# raw cell 1
ls()
savehistory(file="MyHistory.Rhistory") #save all previously run commands in a special formatted file
loadhistory("MyHistory.Rhistory") #load all commands stored in the specified file
hist <- read.delim("MyHistory.Rhistory") #see how the file is formatted: number of line and associated command
head(hist)
## Code cell n°40 ##

x <- c(3, 7, 1, 2) # we define a variable x with 4 numeric values concatenated
x

## Code cell n°41 ##

print(x) 

## Code cell n°42 ##

is.numeric(x)

## Code cell n°43 ##

x < 2 # we test whether the 4 values are < 2

## Code cell n°44 ##

x == 2

## Code cell n°45 ##

class(x)
class(s)
is.character(s)
is.numeric(s)
print(as.numeric(x < 2))
is.numeric("1")
is.numeric(as.numeric("1"))
is.numeric(c(1, "1"))   # c(1, "1") is a concatenation of a numeric value and a character value
class(c(1, "1"))

## Code cell n°46 ##

a <- c()
a

## Code cell n°47 ##

weight <- c(60, 72, 57, 90, 95, 72)
weight

## Code cell n°48 ##

print(weight)

## Code cell n°49 ##

4:10
print(4:10)

## Code cell n°50 ##

print(seq(from = 4, to = 10))
print(seq(4, 10))

## Code cell n°51 ##
? seq

## Code cell n°52 ##

print(seq(from = 2, to = 10, by = 2))

## Code cell n°53 ##

print(rep(x = 4, times = 2))

## Code cell n°54 ##

print(rep(seq(4, 10, 2)))
print(c(rep(1, 4),rep(2, 4)))
print(c(5, s))

## Code cell n°55 ##

class(c(5, s))
length(1:10)
length(weight)
str(weight)

## Code cell n°56 ##

size <- c(1.75, 1.8, 1.65, 1.9, 1.74, 1.69)
print(size^2)
print(bmi <- weight/size^2 )
print(bmi)

## Code cell n°57 ##

print(sort(size))
mean(size)
sd(size)
median(size)
min(size)
max(size)
print(range(size))

summary(size)

## Code cell n°58 ##

print(size)
size[1]
size[2]
size[6]
print(size[c(2,6)])
print(size[c(6,2)])
min(size[c(6,2)])

## Code cell n°59 ##

names(size)
names(size) <- c("Fabien", "Pierre", "Sandrine", "Claire", "Bruno", "Delphine")
size
str(size)

## Code cell n°60 ##

myData <- matrix(c(1, 2, 3, 11, 12, 13), nrow = 2, ncol = 3)
myData
class(myData)

## Code cell n°61 ##

myData <- matrix(c(1, 2, 3, 11, 12, 13), nrow = 2, ncol = 3, byrow = TRUE)
myData

## Code cell n°62 ##

print(dim(myData))
str(myData)
nrow(myData)
ncol(myData)

## Code cell n°63 ##

print(myData)

## Code cell n°64 ##

myData[1, 2] # returns the value of the 1st row and 2nd column

## Code cell n°65 ##

myData[2, 1] # returns the value of the 2nd row and 1st column

## Code cell n°66 ##

print(myData[, 1]) # returns the values of the vector corresponding to the 1st column

## Code cell n°67 ##

print(myData[2,])  # returns the values of the vector corresponding to the 2nd row

## Code cell n°68 ##

myData[, 2:3] # subsets the initial matrix returning a sub-matrix
             # with all rows of the 2nd and 3rd columns from the initial matrix
             # the generated matrix has 2 rows and 2 columns

## Code cell n°69 ##

print(dim(myData[, 2:3])) # the generated matrix has 2 rows and 2 columns

## Code cell n°70 ##

myData[, 1]        # returns the values of the vector corresponding to the 1st column
class(myData[, 1]) # we extract a vector containg numbers -> thus the class is numeric and no more matrix
length(myData[1,])
length(myData[, 1])

## Code cell n°71 ##

myData2 <- cbind(weight, size, bmi) # joining the vectors as columns
myData2
myData3 <- rbind(weight, size, bmi) # joining the vectors as rows
myData3

## Code cell n°72 ##

myData2*2
summary(myData2)
mean(myData2)
mean(myData2[, 1])

## Code cell n°73 ##

myDataf <- data.frame(weight, size, bmi)
myDataf

## Code cell n°74 ##

class(myDataf)

## Code cell n°75 ##

str(myDataf)

## Code cell n°76 ##

print(dim(myDataf))

## Code cell n°77 ##

d <- data.frame()
d
dim(d)

## Code cell n°78 ##

class(myData2)
class(as.data.frame(myData2))
str(as.data.frame(myData2))

## Code cell n°79 ##

d2 <- as.data.frame(cbind(1:2, 10:11))
d2
str(d2)

## Code cell n°80 ##

d <- as.data.frame(matrix(NA, 2, 3))
d
dim(d)
str(d)

## Code cell n°81 ##

path_to_file <- "/shared/projects/2413_rnaseq_cea/alldata/Example_Data/Temperatures.txt" #in case you did not copy it on your personal directory in Pipe06
temperatures <- read.table(path_to_file,
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE)
temperatures
str(temperatures)

## Code cell n°82 ##

temperatures.2 <- read.table(path_to_file,
                             sep = "\t",
                             header = TRUE,
                             stringsAsFactors = TRUE)
str(temperatures.2)

## Code cell n°83 ##

levels(temperatures.2$Month)

## Code cell n°84 ##

# save a dataframe as a text file in the working directory
write.table(myDataf,
            file = "bmi_data.txt",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE)

## Code cell n°85 ##

rm(myDataf)
myDataf <- read.table("bmi_data.txt",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE)
head(myDataf) #myDataf is again accessible
file.remove("bmi_data.txt") #to clean the working directory

## Code cell n°86 ##

print(myDataf)

## Code cell n°87 ##

row.names(myDataf)
names(myDataf)

## Code cell n°88 ##

print(myDataf[, 2])
print(myDataf[, "size"])
print(myDataf$size)

## Code cell n°89 ##

myDataf[2,]

## Code cell n°90 ##

myDataf["Pierre",]

## Code cell n°91 ##

class(myDataf["Pierre",])

## Code cell n°92 ##

temp <- unlist(myDataf["Pierre",])
print(temp)
class(temp)

## Code cell n°93 ##

d2
d2$new <- 1:2
d2

## Code cell n°94 ##

gender <- c("Man", "Man", "Woman", "Woman", "Man", "Woman")
print(gender)
myDataf$sex <- gender
print(myDataf$sex)
myDataf
str(myDataf)

## Code cell n°95 ##

d3 <-  data.frame(d, d2)
d3

## Code cell n°96 ##

d3 <- as.data.frame(rbind(d3, rep("toto", 6)))
d3
str(d3)

## Code cell n°97 ##

print(which ( myDataf$sex == "Woman") )

## Code cell n°98 ##

myDataf [ which ( myDataf$sex == "Woman") , ] 

## Code cell n°99 ##

str(myDataf [ which ( myDataf$sex == "Woman") , ])

## Code cell n°100 ##

print(which ( myDataf$sex != "Man"))

## Code cell n°101 ##
print(which (! myDataf$sex == "Man"))

## Code cell n°102 ##

myDataf2 <- myDataf
myDataf2["Claire", "sex"] <- NA
myDataf2

## Code cell n°103 ##

myDataf2[myDataf2$sex == "Woman",]

## Code cell n°104 ##

myDataf2[which(myDataf2$sex == "Woman"),]

## Code cell n°105 ##

print(grep("Wom", myDataf$sex))

## Code cell n°106 ##

print(grep("Woman", myDataf$sex))

## Code cell n°107 ##

myDataf [grep("Woman", myDataf$sex), ] 

## Code cell n°108 ##

print(grep("a", row.names(myDataf)))

## Code cell n°109 ##

myDataf [grep("a", row.names(myDataf)),]

## Code cell n°110 ##

WomenDataf <- subset(myDataf, gender== "Woman")
WomenDataf

## Code cell n°111 ##

filteredData <- myDataf [ which ( myDataf$sex == "Woman" & myDataf$weight < 80 & myDataf$bmi > 20), ]
filteredData

## Code cell n°112 ##

subset( myDataf, sex == "Woman" & weight < 80 & bmi > 20)

## Code cell n°113 ##

myDataf$index <- 1:6
myDataf

## Code cell n°114 ##

OtherData <- data.frame(c(1:5, 7),rep(c("right-handed","left-handed"),3))
names(OtherData) <- c("ID","handedness")
OtherData

## Code cell n°115 ##

myMergedDataf <- merge(myDataf, OtherData,
                       by.x = "index", by.y = "ID",
                       all.x = TRUE, all.y = TRUE,
                       sort = FALSE)
myMergedDataf

## Code cell n°116 ##

myfolder <- getwd()
#myfolder <- setwd('/shared/ifbstor1/projects/2423_rnaseq_cea/mylogin') # devrait être inutile
#myfolder

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
dir.create(paste0(myfolder,"/Results/Rintro"), recursive = TRUE)

# storing the path to this output folder in a variable
rintrofolder <- paste0(myfolder,"/Results/Rintro/", sep = "")
rintrofolder

save(myDataf, file = paste0(rintrofolder,"myDataf.RData"))

## Code cell 117 ##   

file.copy("/shared/projects/2413_rnaseq_cea/pipeline/Pipe_07b-R_plots_stats_withR.ipynb", myfolder)
## Code cell 118 ##   

# creation of the directory, recursive = TRUE is equivalent to the mkdir -p in Unix
#dir.create(paste0(myfolder,"/run_notebooks"), recursive = TRUE)
runfolder <- paste0(myfolder,"/run_notebooks")
runfolder

file.copy(paste0(myfolder, "/Pipe_07a-R_intro-to-R.ipynb"), paste0(runfolder, "/Pipe_07a-R_intro-to-R_run.ipynb"))