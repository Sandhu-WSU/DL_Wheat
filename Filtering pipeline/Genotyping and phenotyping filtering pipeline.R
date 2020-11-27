# Data filtering pipeline

setwd("")


#install.packages("lme4")
library(lme4)
#install.packages("reshape2")
library("reshape2")
#install.packages("LDcorSV")
library("LDcorSV")
#install.packages("rcompanion")
library(rcompanion)
#install.packages("compiler")
library(compiler)


# Load phenotype data (It contains three environment phenotype data of five agronomic traits for 1944 plants planted at Spillman Farm for 2014-2016. The name of plant is given by taxa and their respective families are also represented.)

myY  <- read.csv("Pheno.csv", header = TRUE)
myY <- myY[order(myY$Taxa, decreasing = FALSE),]
rownames(myY) <- NULL
mynew <- c(1,1:1943)
myY2 <- cbind(myY,mynew)


# Load genetic data which contains 73,345 markers for 2100 Plants

load("myGD")
dim(myGD)


# Load genetic map for 73,345 markers anchored on wheat chromosomes

load("myGM")
dim(myGM)
myGD <- t(myGD)


# Associating phenotype and genotype data

fY  <- myY2[myY2[,2]   %in% rownames(myGD) ,]
fGD <- myGD[rownames(myGD) %in% fY[,2]   ,]
sapply(fY, class)
dat <- transform(fY, Name = factor(Name),
                 Taxa = factor(Taxa),
                 Family = factor(Family),
                 Env   = factor(Env)  ,
                 Yield = as.numeric(Yield),
                 TSTWT = as.numeric(TSTWT),
                 Protein     = as.numeric(Protein),
                 Height = as.numeric(Height),
                 DTH = as.numeric(DTH))
sapply(dat, class)
head(dat)


# Separating the phenotype data of each environment

e2014 <- dat[dat$Env == "2014",c(2,3,6:10)]
head(e2014)
e2015 <- dat[dat$Env == "2015",c(2,6:10)]
head(e2015)
e2016 <- dat[dat$Env == "2016",c(2,6:10)]
head(e2016)
new_dat <- merge(e2014, e2015, by="Taxa")
new_dat <- merge(new_dat, e2016, by="Taxa")
colnames(new_dat) <- c("Taxa","Family", "2014_Yield","2014_TSTWT", "2014_Protein","2014_Height", "2014_DTH", "2015_Yield","2015_TSTWT", "2015_Protein", "2015_Height", "2015_DTH", "2016_Yield", "2016_TSTWT",  "2016_Protein", "2016_Height", "2016_DTH")


# Removing the plants missing phenotyping data

dim(new_dat)
countNAs <- function(X){length(which(is.na(X)))}
new_dat$Missing <- apply(new_dat[,3:17], 1, countNAs)
new_dat$Keep <- ifelse(new_dat$Missing > 0, FALSE, TRUE)
dat_lessNAs <- new_dat[new_dat$Keep == TRUE ,1:17]
rownames(dat_lessNAs) <- NULL
head(dat_lessNAs)
fGD <- myGD[rownames(myGD) %in% dat_lessNAs[,1]   ,]


# Quality check begins
# Remove markers with >20% missing data

M <- fGD
dim(M)
n <- ncol(M)
a <- nrow(M)
TooManyBlanks <- NULL
for (i in 1:n){
  m <- sum(M[,i]==1)
  if (m > a*.2) {TooManyBlanks <- c(TooManyBlanks, i)}
}

length(TooManyBlanks)
dim(TooManyBlanks)
M <- M[,-TooManyBlanks]                         
dim(M)
n <- ncol(M)
a <- nrow(M)


# Check for duplicate markers or individuals

dupIND <- t(M[,sample(seq(1:ncol(M)), 55316, replace=FALSE )])
colnames(dupIND) <- rownames(M)
corIND <- cor(dupIND, use="pairwise.complete.obs")
length(which(corIND==1))  


# Remove markers that are monomorphic 

M[M==1] <- NA       
dim(M)
vars <- matrix(NA, 1, 55316)
for (checkM in 1:55316){                             
  disone <- var(M[,checkM], na.rm = TRUE)
  vars[1,checkM] <- disone
  
}
length(which(vars == 0))                                            
rm_mono <- which(vars == 0)
M <- M[,-rm_mono]
dim(M)
K = dim(rm_mono)
REMOVE_REDUNDANT <- function(M = M_temp){
  
  nMs <- ncol(M)
  corSNP <- cor(M[,1:nMs], use="pairwise.complete.obs")                       
  length(which(corSNP==1))                                                        
  
  index <- which(corSNP==1, arr.ind = TRUE)                                       
  rn <- as.data.frame(rownames(corSNP)[index[,1]])                                
  rn$col <- colnames(corSNP)[index[,2]]                                           
  colnames(rn) <- c("row", "col")
  rn <- transform(rn, row = as.character(row),
                  col = as.character(col))
  rn$count <- apply(rn, 1, function(x){length(unique(x))})                        
  dups <- rn[rn[,3]==2,]                                                          
  dups <- dups[,-3]
  
  Miss_Rate <- matrix(NA, nMs, 1)
  for (mr in 1:nMs) {Miss_Rate[mr,1] <- sum(is.na(M[,mr]))}
  Miss_Rate <- as.data.frame(Miss_Rate)
  Miss_Rate$SNP <- colnames(M)
  
  colnames(dups)[1] <- "SNP"
  theList <- merge(x = dups, y = Miss_Rate, by ="SNP", all.x=TRUE, all.y = FALSE)        
  colnames(theList) <- c("M1", "SNP", "mr1")
  theList <- merge(x = theList, y = Miss_Rate, by ="SNP", all.x=TRUE, all.y = FALSE)    
  colnames(theList) <- c("M2", "M1", "mr1", "mr2")
  theList <- theList[,c(2,1,3,4)]
  
  nr <- nrow(theList)
  Keep <- NULL
  for (thisismylifenow in 1:nr){
    
    if (theList[thisismylifenow,1] %in% theList[(thisismylifenow + 1):nr,1])     
    {next}                                                                       
    #| If M1 exist as a match anywhere else, move to next row
    if (theList[thisismylifenow,1] %in% theList[(thisismylifenow + 1):nr,2])     
    {next}                                                                       
  
    
    if (theList[thisismylifenow,2] %in% theList[(thisismylifenow + 1):nr,1])     
    {next}                                                                       
    #| If M2 exist as a match anywhere else, move to next row
    if (theList[thisismylifenow,2] %in% theList[(thisismylifenow + 1):nr,2])     
    {next}                                                                       
  
    
    if (theList[thisismylifenow,3] > theList[thisismylifenow,4])               
    {Keep <- c(Keep, theList[thisismylifenow,2])}                              
    
    else {Keep <- c(Keep, theList[thisismylifenow,1])}                           
    
  }
  
  Judge <- unique(rbind(as.matrix(theList[,1]), as.matrix(theList[,2])))
  Drop <- Judge[!(Judge %in% Keep)]
  
  return(Drop)
}


# Scan for redundant markers by chromosome and remove those

droppy <- seq(1:21)
Drop_Redundant <- NULL
for (droopy in droppy){
  M_chroms <- myGM[myGM$Chrom == droopy,1]
  M_temp <- M[,which(colnames(M) %in% M_chroms)]
  remove_these <- REMOVE_REDUNDANT(M_temp)
  print(droopy)
  Drop_Redundant <- c(Drop_Redundant ,remove_these)
}
Ms <- colnames(M)
Ms <- Ms[!(Ms %in% Drop_Redundant)]                                             
final <- M[,Ms]                                                                          


# End of >20% Missing/Duplicate Marker Drop     
# Remove plants with >10% missing genotyping data & markers  < 10% MAF

GD_polished <- final
GD_polished[is.na(GD_polished)] <- 1                               
n <- ncol(GD_polished)
a <- nrow(GD_polished)
Missing <- NULL
for (i in 1:a){
  m <- sum(GD_polished[i,]==1)
  if (m > n*.1) {Missing <- c(Missing, i)}
}
dim(missing)
GD_polished <- GD_polished[-Missing,]                             
n <- ncol(GD_polished)
a <- nrow(GD_polished)
TheJudge <- NULL
for (i in 1:n){
  s <- sum(GD_polished[,i]==0)
  if (s < a*0.10 | s > a*0.99) {TheJudge <- c(TheJudge, i)}
}
GD_1MAF <- GD_polished[,-TheJudge]                                  


# Double check new data set for individuals with >10% missing data

n <- ncol(GD_1MAF)
a <- nrow(GD_1MAF)

Missing <- NULL
for (i in 1:a){
  m <- sum(GD_1MAF[i,]==1)
  if (m > n*.1) {Missing <- c(Missing, i)}
}


#| Result: 0 individuals with >20% missing data
# Repeat selection of only common entries

sapply(dat_lessNAs, class)
fY  <- dat_lessNAs[dat_lessNAs[,1]   %in% rownames(GD_1MAF),]
fGD <- GD_1MAF[rownames(GD_1MAF) %in% fY[,1]   ,]
sapply(fY, class)
dat <- fY
rownames(fY) <- fY$Taxa
NAM_dat <- merge(fY, fGD, by="row.names")
GM <- as.matrix(colnames(GD_1MAF))
colnames(GM) <- "SNP"
GM <- merge(x = GM, y = myGM, by = "SNP", all.x = FALSE, all.y = FALSE)
GM <- GM[order(GM$RefSeqv1_position, decreasing = FALSE),]
GM <- GM[order(GM$Chrom, decreasing = FALSE),]
rownames(GM) <- NULL
sapply(NAM_dat[1:10,1:2], factor)
GD <- as.data.frame(NAM_dat[,-c(1,3:18)])
NAM_dat <- NAM_dat[order(NAM_dat$Taxa),]
rownames(NAM_dat) <- NULL


rm(list=setdiff(ls(), c("NAM_dat", "GM", "GD" )))


# This data contains filtered data, GD and GM
save.image("Filtered_Data")


