# Codes for doing genomic selection with ridge regression best linear unbiased predictor model


setwd("")


# This file contains all the filtered phenotypic, genotypic and genetic map data
# Phenotype data is for 635 plants for three environments

load("Filtered_Data")


#install.packages("MASS")
library("MASS")
#install.packages("multtest")
library(multtest)
#install.packages("gplots")
library(gplots)
#install.packages("compiler")
library(compiler)
#install.packages("scatterplot3d")
library("scatterplot3d")
#install.packages("rcompanion")
library(rcompanion)
#install.packages("rrBLUP")
library(rrBLUP)

# Time function will assist in finding how much time each fold takes

timeit <- function(x = TIME, printy = "This segment took") {clock = proc.time() - x
clock = clock [3]
clock = unname(clock)
clock = (clock / 60)
clock = round(clock, digits = 1)

print(printy)
print(paste(clock," minutes"))
}


# Take the phenotype for which you want to make the predictions

Ys <- as.matrix(cbind()) # Enter your phenotypic data here
colnames(Ys) <- c("2016_TSTWT")
head(Ys)


# Select the markers matrix by removing other data columns

M      <- as.matrix(NAM_dat[,-c(1:18)])
head(M[1:5,1:5])


# Random sample number of markers required to be included in the model. In wheat, anything above 4000 gives the same results with rrBLUP
n2 <- ncol(M)
M_sample <- sample(seq(1:n2), 30000, replace=FALSE)
n      <- nrow(Ys)


# Set the required number of replications and folds
nrep  <- 200
nfold <- 5
results_RandomFF <- matrix(NA,nrep*nfold,1)
colnames(results_RandomFF) <- colnames(Ys)


# This function operates on the working principle of rrBLUP package which assumes random marker effect with common variance. Separate training and testing sets are used and GEBVs are predicted for all the individuals.

for (y in 1) {
  
  myY <- as.matrix(Ys[,y])
  
  for (rep in 1:nrep){
    myCut=cut(1:n,nfold,labels=FALSE)      
    fold=sample(myCut,n)                   
    for (f in 1:nfold){   
      
      TIME <- proc.time()
      
      testing  <- (fold==f)                   
      training <- !testing                     
      
      ans.RR <- mixed.solve(y = myY[training,],
                            Z = M[training,M_sample])
      
      ran_effect <- M[,M_sample]     %*% ans.RR$u                       
      GEBV       <-  ran_effect                  
      
      Pred_v_Pheno <- cor(GEBV[testing,],myY[testing,], use="pairwise.complete.obs")
      
      row <- (rep*5) - (5-f)
      results_RandomFF[row,y] <- Pred_v_Pheno
      
      timeit(printy="1 Fold took:")
      
    } 
  }
} 


results_RAN <- as.data.frame(colMeans(results_RandomFF))
RAN <- as.data.frame(results_RandomFF) 
mean(RAN$`2016_TSTWT`)
error_bar<- sd(RAN$`2016_TSTWT`)/sqrt(nrep*nfold)

# We always report the final results as pearson correlation between the observed and predicted values. This will give the idea to plant breeder how efficiently his selections will look depending upon the GEBVs.