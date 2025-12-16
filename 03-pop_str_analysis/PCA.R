#R script PCA
getwd()
rm(list=ls()) ## 
graphics.off() ## 
getwd()
setwd("/share/projects/ANE/denovo_m5M2n3/populations_r75_AZTI_SRA_L1L4_502M/PCA/")
library(adegenet)
library(ggplot2)
library(plyr)
library(dplyr)


data <- read.structure("../plink_filtered_LD02_allSNPs.str", n.ind=419, n.loc=9209, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=1, NA.char=-9, ask=FALSE)
# pop names
levels(data@pop) <- c("ADR", "ALB", "BOB", "CAD", "CAN","ENG_N","ENG_S", "GCA", "IRE", "LYO", "MAR", "POR_N", "POR_S", "UK", "ATL_MAR", "MED_LAG", "ATL_EST", "MED_MAR", "EST", "PLA")
sum(is.na(data$tab))
X <- scaleGen(data, NA.method="mean")
# PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
##  % variance 
PC1 <- (100*pca1$eig/sum(pca1$eig))[1]
PC2 <- (100*pca1$eig/sum(pca1$eig))[2]
PC3 <- (100*pca1$eig/sum(pca1$eig))[3]

PCA_data <- cbind(pca1$li, pop(data))
head(PCA_data)
colnames(PCA_data)<-c("Axis1","Axis2","Axis3","Pop")
write.csv(PCA_data, file="PCA_data_oneSNPperTag_AZTI.csv")

PCA_data_MOD <-read.csv("PCA_data_oneSNPperTag_AZTI_MOD.csv", sep = ";", row.names = 1)
PCA_data_MOD$Pop5 <-factor(PCA_data_MOD$Pop5, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB", "POR_N", "POR_S", "CAD", "MAR", "GCA","ALB", "LIO","ADR")) 
tiff('FigurePCA_9209SNPs_419indv_oneSNPperTAG_PC1_PC2.tiff', width = 20, height = 15, units = "cm", res=400)
ggplot(data = PCA_data_MOD, aes(x=Axis1, y=Axis2, color=Pop5))+ geom_point(aes(shape=Shape), fill="#fafafa",stroke = 1,  size = 3) + theme_bw()+ stat_ellipse(level=0.75, linetype = "dotted") +
  scale_color_manual(values=c("#8ff7de", "#1fb487", "#216444", "yellowgreen", "#4775ff","#8b0023", "#ff0000","#ff8181", "#d87a00", "gold", "#bc92ea", "#8600b3","#fc0fc0"))
  scale_shape_manual(values=c(21, 19)) + theme(legend.title= element_blank())+  theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = 12,face = "bold"))+
labs(title = "9209SNPs_419indv", x ="PC1 (X%)", y = "PC2 (X%)") 
dev.off()
