# title: "PCA scores and PLS"

## packages
library(dplyr)
library(mixOmics)
library(tidyr)
library(ggraph)
library(igraph)
library(tidyverse)
library(MUVR)
library(doParallel)

## Data loading
sperm <- read.csv("dataset_journal.csv", check.names = FALSE)
metadata <- sperm[,2:16]
colnames(metadata) <- sperm[,1]
metabolomics <- sperm[,17:length(sperm)]


## Normalizacion metabolomics
normalize<-function(data,vars=colnames(data),group_var,method){
  subdata<-data[,vars]
  subdata<-switch(method,
                  auto=apply(subdata,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)),        
                  level=apply(subdata,2,function(x) (x-mean(x,na.rm=T))/mean(x,na.rm=T)),
                  log=apply(subdata,2,function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sd(log10(x+1),na.rm=T)),
                  vast=apply(subdata,2,function(x) ((x-mean(x,na.rm=T))/sd(x,na.rm=T))*(mean(x,na.rm=T)/sd(x,na.rm=T))),         
                  logpareto=apply(subdata,2,function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sqrt(sd(log10(x+1),na.rm=T))),
                  pareto=apply(subdata,2,function(x) (x-mean(x,na.rm=T))/sqrt(sd(x),na.rm=T))
  )
  data[,vars]<-subdata
  return(data)
}
metabolomicsNorm <- normalize(metabolomics, method = "logpareto")

#####################
# PCA score vectors #
#####################

## Sperm Quality vector
quality <- c("Viable sperm", "Progressive Motility", "Total Motility", "Sperm w/ normal morphology")
meta_quality <- metadata[, quality]
pca_quality <- mixOmics::pca((meta_quality), ncomp = 2, center = TRUE, scale = TRUE)
pca_quality$loadings # For supplementary Table 1

## Sperm Function vector
function. <- c("Intracellular calcium levels", "Sperm mitochondrial potential", "Viable sperm w/ intact acrosome")
meta_function <- metadata[, function.]
meta_function <- meta_function[meta_function$Id != "FIV 3 10",] ## Removing sample with outlier




pca_function <- mixOmics::pca((meta_function), ncomp = 1, center = TRUE, scale = TRUE)
pca_function$loadings

## IVF outcomes vector
IVF <- c("Early blastocyst/blastocyst", "Total embryos", "Hatching/Hatched blastocyst", "Morulae", "Morulae and total blastocyst", "Developmental competency fertilised oocyte", "Developmental ratio", "Fertilisation Rate")
meta_IVF <- metadata[, IVF]
pca_IVF <- mixOmics::pca((meta_IVF), ncomp = 2, center = TRUE, scale = TRUE)
pca_IVF$loadings

#################
# PLS modelling #
#################

## PLS Sperm Quality
scores.quality <- pca_quality[["variates"]][["X"]]
cl <- parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
model_quality <- MUVR::MUVR(X = (metabolomicsNorm), 
                            Y = scores.quality[,"PC1"], 
                            ML = F, 
                            method ='PLS',
                            nRep = 30, 
                            nOuter = 4, 
                            varRatio = 0.9,
                            scale = FALSE)
stopCluster(cl)
plotMV(model_quality) ## For Figure 1
sig_LR_PLS_Quality <- MUVR::getVIP(model_quality, model = "max")

## PLS Sperm Function
metabolomicsNormReduced <- metabolomicsNorm[metabolomicsNorm$Id != "FIV 3 10",] ## Removing sample with outlier
scores.function <- pca_function[["variates"]][["X"]]
cl <- parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
model_function <- MUVR::MUVR(X = (metabolomicsNormReduced), 
                            Y = scores.function[,"PC1"], 
                            ML = F, 
                            method ='PLS',
                            nRep = 30, 
                            nOuter = 4, 
                            varRatio = 0.9,
                            scale = FALSE)
stopCluster(cl)
plotMV(model_function) ## For Figure 1
sig_LR_PLS_function <- MUVR::getVIP(model_function, model = "max")

## PLS IVF
scores.IVF <- pca_IVF[["variates"]][["X"]]
cl <- parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
model_IVF <- MUVR::MUVR(X = (metabolomicsNorm), 
                            Y = scores.IVF[,"PC1"], 
                            ML = F, 
                            method ='PLS',
                            nRep = 30, 
                            nOuter = 4, 
                            varRatio = 0.9,
                            scale = FALSE)
stopCluster(cl)
plotMV(model_IVF) ## For Figure 1
sig_LR_PLS_IVF <- MUVR::getVIP(model_IVF, model = "max")

####################
# Permutation test #
####################

cl = parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
permQuality <- MUVR::permutations(model_quality, nPerm = 500)
MUVR::pPerm(model_quality[["fitMetric"]][["Q2"]][2], permQuality)
MUVR::plotPerm(model_quality[["fitMetric"]][["Q2"]][2], permQuality, type ="t",xlab = "Q2",side = "greater", breaks = 50, pos = 3)
stopCluster(cl)

cl = parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
permFunction <- MUVR::permutations(model_Function, nPerm = 500)
MUVR::pPerm(model_Function[["fitMetric"]][["Q2"]][2], permFunction)
MUVR::plotPerm(model_Function[["fitMetric"]][["Q2"]][2], permFunction, type ="t",xlab = "Q2",side = "greater", breaks = 50, pos = 3)
stopCluster(cl)

cl = parallel::makeCluster(parallel::detectCores()-1)       
doParallel::registerDoParallel(cl)
permIVF <- MUVR::permutations(model_IVF, nPerm = 500)
MUVR::pPerm(model_IVF[["fitMetric"]][["Q2"]][2], permIVF)
MUVR::plotPerm(model_IVF[["fitMetric"]][["Q2"]][2], permIVF, type ="t",xlab = "Q2",side = "greater", breaks = 50, pos = 3)
stopCluster(cl)

