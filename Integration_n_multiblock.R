## Packages
library(mixOmics)
library(dplyr)
library(tidyr)

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

## N-block_PLS
quality <- c("Viable sperm", "Progressive Motility", "Total Motility", "Sperm w/ normal morphology")
meta_quality <- metadata[, quality]
function. <- c("Intracellular calcium levels", "Sperm mitochondrial potential", "Viable sperm w/ intact acrosome")
meta_function <- metadata[, function.]
IVF <- c("Early blastocyst/blastocyst", "Total embryos", "Hatching/Hatched blastocyst", "Morulae", "Morulae and total blastocyst", "Developmental competency fertilised oocyte", "Developmental ratio", "Fertilisation Rate")
meta_IVF <- metadata[, IVF]
X <- list(Metabolomics = as.matrix(metabolomicsNorm),
          `Sperm Quality` = as.matrix(meta_quality),
          `Sperm Function` = as.matrix(meta_function),
          `Fertilizing Capacity` = as.matrix(meta_IVF)
)

block.pls.result <- block.pls(X, indY = 4, scale = TRUE, mode = "regression")
plotIndiv(block.pls.result, ind.names = FALSE, legend=TRUE, title = "")

## Components correlation for supplementary Table 2
cor1 <- cor(block.pls.result$variates$Metabolomics, block.pls.result$variates$`Sperm Quality`)
rownames(cor1) <- c("comp1_metabolomics","comp2_metabolomics")
colnames(cor1) <- c("comp1_Sperm_Quality","comp2_Sperm_Quality")

cor2 <- cor(block.pls.result$variates$Metabolomics, block.pls.result$variates$`Sperm Function`) 
rownames(cor2) <- c("comp1_metabolomics","comp2_metabolomics")
colnames(cor2) <- c("comp1_Sperm_Function","comp2_Sperm_Function")

cor3 <- cor(block.pls.result$variates$Metabolomics, block.pls.result$variates$`Fertilizing Capacity`)
rownames(cor3) <- c("comp1_metabolomics","comp2_metabolomics")
colnames(cor3) <- c("comp1_Fertilizing_Capacity","comp2_Fertilizing_Capacity")

cor4 <- cor(block.pls.result$variates$`Sperm Quality`, block.pls.result$variates$`Sperm Function`)
rownames(cor4) <- c("comp1_Sperm_Quality","comp2_Sperm_Quality")
colnames(cor4) <- c("comp1_Sperm_Function","comp2_Sperm_Function")

cor5 <- cor(block.pls.result$variates$`Sperm Quality`, block.pls.result$variates$`Fertilizing Capacity`)
rownames(cor5) <- c("comp1_Sperm_Quality","comp2_Sperm_Quality")
colnames(cor5) <- c("comp1_Fertilizing_Capacity","comp2_Fertilizing_Capacity")

cor6 <- cor(block.pls.result$variates$`Sperm Function`, block.pls.result$variates$`Fertilizing Capacity`)
rownames(cor6) <- c("comp1_Sperm_Function","comp2_Sperm_Function")
colnames(cor6) <- c("comp1_Fertilizing_Capacity","comp2_Fertilizing_Capacity")


## Creation of the network
corPlot <- network(block.pls.result, blocks = c(1,2,3,4), cutoff = 0.3, save = 'jpeg', name.save = 'blocPLS')

## Extracting correlation matrix from the network
M_Metabolomics_Sperm_Quality <- as.data.frame(corPlot$`M_Metabolomics_Sperm Quality`)
M_Metabolomics_Sperm_Function <- as.data.frame(corPlot$`M_Metabolomics_Sperm Function`)
M_Metabolomics_Fertilizing_Capacity <- as.data.frame(corPlot$`M_Metabolomics_Fertilizing Capacity`)
M_Sperm_Quality_Sperm_Function <- as.data.frame(corPlot$`M_Sperm Quality_Sperm Function`)
M_Sperm_Quality_Fertilizing_Capacity <- as.data.frame(corPlot$`M_Sperm Quality_Fertilizing Capacity`)
M_Sperm_Function_Fertilizing_Capacity <- as.data.frame(corPlot$`M_Sperm Function_Fertilizing Capacity`)

M_Metabolomics_Sperm_Quality$from <- rownames(M_Metabolomics_Sperm_Quality)
M_Metabolomics_Sperm_Function$from <- rownames(M_Metabolomics_Sperm_Function)
M_Metabolomics_Fertilizing_Capacity$from <- rownames(M_Metabolomics_Fertilizing_Capacity)
M_Sperm_Quality_Sperm_Function$from <- rownames(M_Sperm_Quality_Sperm_Function)
M_Sperm_Quality_Fertilizing_Capacity$from <- rownames(M_Sperm_Quality_Fertilizing_Capacity)
M_Sperm_Function_Fertilizing_Capacity$from <- rownames(M_Sperm_Function_Fertilizing_Capacity)

## Gathering matrizes into long column. Supplementary Table 3
connect_1 <- M_Metabolomics_Sperm_Quality %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_1 <- connect_1[connect_1$value != 0,]
connect_2 <- M_Metabolomics_Sperm_Function %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_2 <- connect_2[connect_2$value != 0,]
connect_3 <- M_Metabolomics_Fertilizing_Capacity %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_3 <- connect_3[connect_3$value != 0,]
connect_4 <- M_Sperm_Quality_Sperm_Function %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_4 <- connect_4[connect_4$value != 0,]
connect_5 <- M_Sperm_Quality_Fertilizing_Capacity %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_5 <- connect_5[connect_5$value != 0,]
connect_6 <- M_Sperm_Function_Fertilizing_Capacity %>% pivot_longer(!from, names_to = "to", values_to = "value")
connect_6 <- connect_6[connect_6$value != 0,]

connect_Z <- rbind(connect_1, connect_2, connect_3, connect_4, connect_5, connect_6)
