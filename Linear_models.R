# title: "Linear models"

## packages


## Data loading
sperm <- read.csv("dataset_journal.csv", check.names = FALSE)
metadata <- sperm[,2:16]
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

## Imputing META (zeros)
metadata[metadata == 0] <- NA

for(i in 1:ncol(metadata)){
  metadata[is.na(metadata[,i]), i] <- min(metadata[,i], na.rm = TRUE)/2}

#####################
# PCA score vectors #
#####################

## Sperm Quality vector
quality <- c("Viable sperm", "Progressive Motility", "Total Motility", "Sperm w/ normal morphology")
meta_quality <- metadata[, quality]
pca_quality <- mixOmics::pca(log(meta_quality), ncomp = 2, center = TRUE, scale = TRUE)
pca_quality$loadings # For supplementary Table 1

## Sperm Function vector
function. <- c("Intracellular calcium levels", "Sperm mitochondrial potential", "Viable sperm w/ intact acrosome")
meta_function <- metadata[, function.]
pca_function <- mixOmics::pca(log(meta_function), ncomp = 2, center = TRUE, scale = TRUE)
pca_function$loadings

## IVF outcomes vector
IVF <- c("Early blastocyst/blastocyst", "Total embryos", "Hatching/Hatched blastocyst", "Morulae", "Morulae and total blastocyst", "Developmental competency fertilised oocyte", "Developmental ratio", "Fertilisation Rate")
meta_IVF <- metadata[, IVF]
pca_IVF <- mixOmics::pca(log(meta_IVF), ncomp = 2, center = TRUE, scale = TRUE)
pca_IVF$loadings

#################
# Linear models #
#################

## Sperm Quality
pca_quality <- mixOmics::pca(log(meta_quality), ncomp = 2, center = FALSE, scale = TRUE)
scores.quality <- pca_quality[["variates"]][["X"]]
t <- scores.quality[,"PC1"]
stat <- function(x){(lm(t ~ x, metabolomics))}
abcd <- apply(data.matrix(metabolomics), 2, stat)
p_quality <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Pr(>|t|)"]}))
names(p_quality) <- colnames(metabolomicsNorm)
p_quality <- as.data.frame(p_quality)
p_quality$name = row.names(p_quality)
result <- p_quality
fdr.quality <- p.adjust(result$p_quality, method = "fdr")
result <- cbind(result, fdr.quality)
slope_quality <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Estimate"]}))
names(slope_quality) <- colnames(metabolomicsNorm)
slope_quality <- as.data.frame(slope_quality)
slope_quality$name = row.names(slope_quality)
result <- merge(result, slope_quality, by = "name")

## Sperm Function
pca_function <- mixOmics::pca(log(meta_function), ncomp = 2, center = TRUE, scale = TRUE)
scores.function <- pca_function[["variates"]][["X"]]
t <- scores.function[,"PC1"]
stat <- function(x){(lm(t ~ x, metabolomics))}
abcd <- apply(data.matrix(log(metabolomics)), 2, stat)
p_function <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Pr(>|t|)"]}))
names(p_function) <- colnames(metabolomicsNorm)
p_function <- as.data.frame(p_function)
p_function$name = row.names(p_function)
result <- merge(result, p_function, by = "name")
fdr.calidad <- p.adjust(result$p_function, method = "fdr")
result <- cbind(result, fdr.calidad)
slope_function <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Estimate"]}))
names(slope_function) <- colnames(metabolomicsNorm)
slope_function <- as.data.frame(slope_function)
slope_function$name = row.names(slope_function)
result <- merge(result, slope_function, by = "name")

## IVF
pca_IVF <- mixOmics::pca(log(meta_IVF), ncomp = 2, center = FALSE, scale = TRUE)
scores.IVF <- pca_IVF[["variates"]][["X"]]
t <- scores.IVF[,"PC1"]
stat <- function(x){(lm(t ~ x, metabolomics))}
abcd <- apply(data.matrix(log(metabolomics)), 2, stat)
p_IVF <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Pr(>|t|)"]}))
names(p_IVF) <- colnames(metabolomicsNorm)
p_IVF <- as.data.frame(p_IVF)
p_IVF$name = row.names(p_IVF)
result <- merge(result, p_IVF, by = "name")
fdr.IVF <- p.adjust(result$p_IVF, method = "fdr")
result <- cbind(result, fdr.IVF)
slope_IVF <- unlist(lapply(abcd, function(x) {summary(x)[["coefficients"]]["x","Estimate"]}))
names(slope_IVF) <- colnames(metabolomicsNorm)
slope_IVF <- as.data.frame(slope_IVF)
slope_IVF$name = row.names(slope_IVF)
result <- merge(result, slope_IVF, by = "name")