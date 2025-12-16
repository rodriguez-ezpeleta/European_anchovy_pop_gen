
##PCAadapt
rm(list=ls()) ## 
graphics.off() ## 
setwd("/share/projects/ANE/denovo_m5M2n3/populations_r75_AZTI_SRA_L1L4_502M/pcadapt")
library(pcadapt)

### READING DATA ###

filename <- read.pcadapt("plink_filtered_LD02_allSNPs.ped", type = "ped", type.out = "matrix") 
popmap <- read.table("../plink_filtered_LD02_allSNPs.nosex", stringsAsFactors = TRUE)
names(popmap) <- c("Pop","Samples")
str(popmap)
popmap$Pop <- factor(popmap$Pop, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB", "CAN", "POR_N", "POR_S", "CAD", "MAR", "GCA", "ALB", "LYO", "ADR", "PLA","EST", "ATL_EST", "MED_LAG", "ATL_MAR", "MED_MAR"))
col <- c("#8ff7de", "#1fb487", "#216444", "yellowgreen", "#4775ff", "#4775ff", "#8b0023", "#ff0000","#ff8181", "#d87a00", "gold", "#bc92ea", "#8600b3","#fc0fc0", "#4775ff", "#4775ff", "#4775ff", "#8600b3", "#4775ff", "#8600b3")

### CHOOSING THE NUMBER OF K OF PRINCIPAL COMPONENTS ###
x <- pcadapt(input = filename, K = 20)
#tiff("PCAdapt_screeplot_LD02.tiff", units="in", width=5, height=5, res=300)
plot(x, option = "screeplot")
#dev.off()
plot(x, option = "scores", pop = popmap$Pop, col=col)
plot(x, option = "scores", i = 1, j = 3, pop = popmap$Pop, col=col)
#The third and the fourth principal components do not ascertain population structure anymore.
plot(x, option = "scores", i = 2, j = 3, pop = popmap$Pop, col=col)## No more structure
plot(x, option = "scores", i = 3, j = 4, pop = popmap$Pop, col=col)## No more structure
#Use the screeplot to define how many PC you need to discover outlier SNPs
x <- pcadapt(filename, K = 3)
summary(x)


### GRAPHICAL TOOLS ###

## Manhattan Plot
#tiff("PCAdapt_manhattanplot_K3_LD02.tiff", units="in", width=5, height=5, res=300)
plot(x, option = "manhattan")
#dev.off()
## Q-Q Plot
  #Check the expected uniform distribution of the p-values using a Q-Q plot. 
plot(x, option = "qqplot")
## Histograms of the test statistic and of the p-values
    #An histogram of p-values confirms that most of the p-values follow an uniform distribution. The excess of small p-values indicates the presence of outliers.
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
    #The presence of outliers is also visible when plotting a histogram of the test statistic Dj.
plot(x, option = "stat.distribution")

### CHOOSING A CUTOFF FOR OUTLIER DETECTION ###
#q-values: The R package qvalue, transforms p-values into q-values. For a given alfa (real valued number between 0 and 1), SNPs with q-values less than alfa will be considered as outliers with an expected false discovery rate bounded by alfa. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. For an expected false discovery rate lower than 5%:
#Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
write.table(outliers,"Numbers_Benjamini_Hochberg_10_LD02.txt", sep="\t")
list_snps <- read.table("plink_filtered_LD02_allSNPs.map")[,1:2]
outliers_snps <- list_snps[outliers,]
neutral_snps <- list_snps[-outliers,]

write.table(outliers_snps[,2], "pcadapt_Benjamini_outlier_snps_AZTI_SRA_L1L4_502M_anchovy_LD02.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(neutral_snps[,2], "pcadapt_Benjamini_neutral_snps_AZTI_SRA_L1L4_502M_anchovy_LD02.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


### PLOTS ###

#1. PCA WITH SNPs UNDER SELECTION #

data <- read.structure("outlier_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.str", n.ind=419, n.loc=716, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=1, NA.char=-9, ask=FALSE)
# pop names
levels(data@pop) <- c("ADR", "ALB", "BOB", "CAD", "CAN","ENG_N","ENG_S", "GCA", "IRE", "LYO", "MAR", "POR_N", "POR_S", "UK", "ATL_MAR", "MED_LAG", "ATL_EST", "MED_MAR", "EST", "PLA")
sum(is.na(data$tab))
X <- scaleGen(data, NA.method="mean")
# PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
##  % of variance
PC1 <- (100*pca1$eig/sum(pca1$eig))[1]
PC2 <- (100*pca1$eig/sum(pca1$eig))[2]
PC3 <- (100*pca1$eig/sum(pca1$eig))[3]
PCA_data <- cbind(pca1$li, pop(data))
head(PCA_data)
colnames(PCA_data)<-c("Axis1","Axis2","Axis3","Pop")
write.csv(PCA_data, file="Outlier_PCA_data_oneSNPperTag_AZTI_LD02.csv")

# I have changed the pop labels in excel
PCA_data_MOD <-read.csv("Outlier_PCA_data_oneSNPperTag_AZTI_LD02_MOD.csv", sep = ";", row.names = 1)
PCA_data_MOD$Pop5 <-factor(PCA_data_MOD$Pop5, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB", "POR_N", "POR_S", "CAD", "MAR", "GCA","ALB", "LIO","ADR")) 
tiff('outlier_pcadapt_716_SNPs_419indv_oneSNPperTAG_LD02_PC1_PC2.tiff', width = 20, height = 15, units = "cm", res=400)
ggplot(data = PCA_data_MOD, aes(x=Axis1, y=Axis2, color=Pop5))+ geom_point(aes(shape=Shape), fill="#fafafa",stroke = 1,  size = 3) + theme_bw()+ stat_ellipse(level=0.75, linetype = "dotted") +
  scale_color_manual(values=c("#8ff7de", "#1fb487", "#216444", "yellowgreen", "#4775ff","#8b0023", "#ff0000","#ff8181", "#d87a00", "gold", "#bc92ea", "#8600b3","#fc0fc0")) +
  scale_shape_manual(values=c(21, 19)) + theme(legend.title= element_blank())+ theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = 12,face = "bold"))+
labs(title = "716SNPs 419indv", x ="PC1 (X%)", y = "PC2 (X%)")
dev.off()

#### PCA WITH NEUTRAL SNPs ####

data <- read.structure("neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.str", n.ind=419, n.loc=6498, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=1, NA.char=-9, ask=FALSE)
# pop names
levels(data@pop) <- c("ADR", "ALB", "BOB", "CAD", "CAN","ENG_N","ENG_S", "GCA", "IRE", "LYO", "MAR", "POR_N", "POR_S", "UK", "ATL_MAR", "MED_LAG", "ATL_EST", "MED_MAR", "EST", "PLA")
sum(is.na(data$tab))
X <- scaleGen(data, NA.method="mean")
# PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
##  % of variance
PC1 <- (100*pca1$eig/sum(pca1$eig))[1]
PC2 <- (100*pca1$eig/sum(pca1$eig))[2]
PC3 <- (100*pca1$eig/sum(pca1$eig))[3]
PCA_data <- cbind(pca1$li, pop(data))
head(PCA_data)
colnames(PCA_data)<-c("Axis1","Axis2","Axis3","Pop")
write.csv(PCA_data, file="Neutral_PCA_data_oneSNPperTag_AZTI_LD02.csv")

# I have changed the pop labels in excel
PCA_data_MOD <-read.csv("Neutral_PCA_data_oneSNPperTag_AZTI_LD02_MOD.csv", sep = ";", row.names = 1)
PCA_data_MOD$Pop5 <-factor(PCA_data_MOD$Pop5, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB", "POR_N", "POR_S", "CAD", "MAR", "GCA","ALB", "LIO","ADR")) 
tiff('neutral_pcadapt_6498SNPs_419indv_oneSNPperTAG_LD02_PC1_PC2.tiff', width = 20, height = 15, units = "cm", res=400)
ggplot(data = PCA_data_MOD, aes(x=Axis1, y=Axis2, color=Pop5))+ geom_point(aes(shape=Shape), fill="#fafafa",stroke = 1,  size = 3) + theme_bw()+ stat_ellipse(level=0.75, linetype = "dotted") +
  scale_color_manual(values=c("#8ff7de", "#1fb487", "#216444", "yellowgreen", "#4775ff","#8b0023", "#ff0000","#ff8181", "#d87a00", "gold", "#bc92ea", "#8600b3","#fc0fc0")) +
  scale_shape_manual(values=c(21, 19)) + theme(legend.title= element_blank())+  theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = 12,face = "bold"))+
  labs(title = "6498SNP 419indv", x ="PC1 (X%)", y = "PC2 (X%)")
dev.off()
