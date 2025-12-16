
#### in R #####
## After running "ADMIXTURE_01.sh"
getwd()
rm(list=ls()) 
graphics.off() 
setwd("/share/projects/ANE/denovo_m5M2n3/populations_r75_AZTI_SRA_L1L4_502M/ADMIXTURE")
getwd()
list.files()
## 
library(vcfR)
library(adegenet)
library(ggplot2)
library(readxl)
library(ape)
library(dplyr)
library(wrapr)
library(tidyverse)
library(tidyr)

#plink19 --file neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02 --recode vcf --out neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02
#plink19 --file outlier_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02 --recode vcf --out outlier_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02
#cat neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.ped | cut -d ' ' -f2 > samples.txt
#cat neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.ped | cut -d ' ' -f1 > pops.txt
#paste samples.txt pops.txt > meta.txt
#rm samples.txt pops.txt

## LOAD DATA
anchovy_vcf <- read.vcfR("neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.vcf")
# create genind
anchovy_data <- vcfR2genind(anchovy_vcf)
# import pop info and other
infopop <- read.table("meta.txt", header = FALSE, sep = "\t")
infopop$x <- paste(infopop$V2,infopop$V1, sep = "_")
#check matching samples names
all(rownames(anchovy_data@tab) == infopop$x)
anchovy_data@pop <- as.factor(infopop$V2)
rownames(anchovy_data@tab) <- as.factor(infopop$V1)

## CHECK K DISTRIBUTION
log<-read.table("./cross_validation/cross_validation_neutral.txt")[,c(3:4)]
log$V3<-gsub("\\(K=", "", log$V3)
log$V3<-gsub("):", "", log$V3)
log$V3<-as.numeric(log$V3)
colnames(log)<-c("Kvalue","cross.validation.error")
#make plot showing the cross validation error across K values 1:10
tiff("ADX_cross_validation_plot_neutral_LD.tiff", units="in", width=10, height=5, res=300)
ggplot(data=log, aes(x=Kvalue, y=cross.validation.error, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point(size = 3) +
  ylab("cross-validation error")+
  xlab("K")+
  theme_classic()+
  scale_x_continuous(breaks = c(1:20))+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))
dev.off()


## PREPARE DATA FOR ADMIXTURE PLOT
#K=3 (do this for each K)
admix<-read.table("./cross_validation/neutral_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02.3.Q")
Pop_ID <- anchovy_data@pop
Samples_names <- rownames(anchovy_data@tab)
admix3 <- cbind(Samples_names,Pop_ID,admix)
colnames(admix3)<-c("Sample","Pop","K1","K2","k3")
admix3_table <-gather(admix3, key="K", value="Prob", K1:K3)
write.csv(admix3_table,"vcf_admix3.csv")
# I have modified the pop labels in excel
admix3_table$Pop <- factor(admix3_table$Pop, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB","BOB_EST", "POR_N", "POR_S", "CAD", "MAR", "GCA", "ALB", "LIO", "LIO_EST", "ADR"))
admix3_table_MOD <- admix3_table %>% mutate(across('Pop', str_replace, 'BOB_EST', 'BOB.'))
admix3_table_MOD <- admix3_table_MOD %>% mutate(across('Pop', str_replace, 'LIO_EST', 'LIO.'))
admix3_table_MOD$Pop <- factor(admix3_table_MOD$Pop, levels=c("IRE", "UK", "ENG_N", "ENG_S", "BOB","BOB.", "POR_N", "POR_S", "CAD", "MAR", "GCA", "ALB", "LIO", "LIO.", "ADR"))

# (optional) order data for the plot:
ad3<-admix3_table[order(admix3_table_MOD$Pop, admix3_table$K, as.numeric(admix3_table$Prob)),]



## PLOT
tiff("ADX_k3_neutral.tiff", units="in", width=15, height=5, res=300)
ggplot(data=admix3_table_MOD,aes(y=Prob, x=fct_inorder(Sample), fill = factor(K)))+
geom_bar(position="fill", stat="identity", width = 1.2)+
scale_fill_manual(values=c("#1e1e1e", "#cecece", "#767676"))+
facet_grid(~Pop, scales = "free",space="free_x") +
theme(axis.title.y = element_blank(),
        axis.text.x =  element_blank(),
        axis.text.y = element_text(size=10, face = "bold"), 
        axis.title.x = element_blank(),
        strip.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold", size=10, angle = 90))+
        theme(legend.position = "none")+
  theme(panel.spacing = unit(0.1, "lines"))
dev.off()


