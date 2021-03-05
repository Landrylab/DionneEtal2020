source("http://bioconductor.org/biocLite.R")
biocLite("EBImage")
install.packages('gitter')
require(gitter)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(magrittr)
library(stringr)
library(UpSetR)
library(gitter)
library(gplots)
library(RColorBrewer)

# Script for the analysis of the second PCA screen, including Sla1 SH3-shuffled baits.

# Gitter Analysis of S2 plate pictures.

Ref="C:/Users/ugols/Desktop/PCA_Fev2019/S2_Pictures_Fev19/IMG_03249.JPG"
f="C:/Users/ugols/Desktop/PCA_Fev2019/S2_Pictures_Fev19"
path_dat="C:/Users/ugols/Desktop/PCA_Fev2019/S2_Pictures_Fev19"
path_grid="C:/Users/ugols/Desktop/PCA_Fev2019/S2_Pictures_Fev19"
gitter.batch(image.files=f,ref.image.file=Ref,plate.format=1536, verbose='p',grid.save=path_grid, dat.save=path_dat)

# Gitter Analysis of DHFR-PCA plate pictures.

Ref="C:/Users/ugols/Desktop/PCA_Fev2019/MTX2_Last_Day_Fev19/d000103_300_001_19-02-13_17-51-46.JPG"
f="C:/Users/ugols/Desktop/PCA_Fev2019/MTX2_Last_Day_Fev19"
path_dat="C:/Users/ugols/Desktop/PCA_Fev2019/Gitter_Analysis_S2_MTX2Last_Fev19"
path_grid="C:/Users/ugols/Desktop/PCA_Fev2019/Gitter_Analysis_S2_MTX2Last_Fev19"
gitter.batch(image.files=f,ref.image.file=Ref,plate.format=1536, verbose='p',grid.save=path_grid, dat.save=path_dat)

# Merge Plate definitions, Prey positions on array with the Gitter analysis.

platesDHFR12_def_Fev2019 = read.csv("Plate_Baits_PCAFev2019_UD.csv", header = T, sep = ";")
platesDHFR3_def_Fev2019 = read.csv("PreysArray_Fev2019.csv", header = T, sep = ",")

PCA19 = data.frame(matrix(vector(), 0, 10, dimnames = list(c(), c("row", "col", "size", "circularity", "flags", "plate", "DHFR1.2", "condition", "rep", "DHFR3"))))
head(PCA19)
for (i in 1:nrow(platesDHFR12_def_Fev2019)) {
  data = read.table(paste("Gitter_Analysis_S2_MTX2Last_Fev19/",as.character(platesDHFR12_def_Fev2019[i,"dat_Gitter"]), sep=""), sep="\t", header=F, stringsAsFactor=F, strip.white = T,comment.char = "#", na.strings="")
  names(data) <- c("row", "col", "size", "circularity", "flags")
  plate = platesDHFR12_def_Fev2019[i, "plate"]
  DHFR1.2 = platesDHFR12_def_Fev2019[i, "DHFR1.2"]
  condition = platesDHFR12_def_Fev2019[i, "condition"]
  rep = platesDHFR12_def_Fev2019[i, "rep"]
  subdf <- cbind(data, rep(plate, nrow(data)), rep(DHFR1.2, nrow(data)), rep(condition, nrow(data)), rep(rep, nrow(data)))
  subdf <- merge(subdf, platesDHFR3_def_Fev2019[,c("col", "row", "orf")], by.x=c("col", "row"), by.y=c("col", "row"))
  names(subdf) <- c("row", "col", "size", "circularity", "flags", "plate", "DHFR1.2", "condition", "rep", "DHFR3")
  PCA19 = rbind(PCA19, subdf)
}
write.csv(PCA19,"AlldataFev2019.csv")

# The dataset available with the article of the raw data for the second DHFR-PCA screen correspond to this step. Some colunm names and some baits names have been edited manually for the paper.

PCA19$size = as.numeric(as.character(PCA19$size))
PCA19$col = as.numeric(as.character(PCA19$col))
PCA19$row = as.numeric(as.character(PCA19$row))
PCA19$plate = as.numeric(as.character(PCA19$plate))
PCA19$rep = as.numeric(as.character(PCA19$rep))
PCA19$DHFR3 = as.character(PCA19$DHFR3)
PCA19$DHFR1.2 = as.character(PCA19$DHFR1.2)

# Analysis of the colonies from the diploid selection S2. Remove colonies that didn't grow after the diploid selection and then remove the S2 data.

S2= subset(PCA19, PCA19$condition == "S2")
S2null = subset(PCA19, PCA19$condition == "S2" & PCA19$size < 500)
PCA19 = anti_join(PCA19, S2null, by = c("row", "col", "plate"))
PCA19 %<>% filter(condition!="S2")

# Log2 transformation of colony size and remove positions with no growth in the DHFR-PCA selection.

PCA19 %<>% mutate(adjust = size+1)
PCA19 %<>% mutate(log2 = log2(adjust))
PCA19 %<>% mutate(baitXreplicate = interaction(DHFR1.2,rep))
PCA19 = subset(PCA19, PCA19$log2 > 2)

# Normalize the colony size to plate background.

PCA19 %<>% 
  group_by(baitXreplicate) %>%
  mutate(med = median(log2, na.rm=T)) %>%
  ungroup()
PCA19 %<>% mutate(mednorm = log2 - med)

# Remove borders.

PCA19SLA1%<>% filter(DHFR3!="border") 

# Distribution of plate background.

disMed = ggplot(data=PCA19SLA1, aes(x = med, colour = DHFR1.2)) +geom_density() + geom_vline(xintercept=0)
pdf("distributionMed_1219Analysis.pdf",7,10)
disMed + theme(legend.position="none")
dev.off()

# 3 plates have almost no background level which results in high plate median (colonies with no growth were removed). They are the only plates with plate median above 8.
# The 3 plates are the second plate replicate of the three WT re-insertion controls. They are assigned the plate median corresponding to the average of first plate replicate of these three baits (6.20, 6.50, 6.35). The assigend plate median is 6.35.

PCA19SLA1LowBKMedOK = subset(PCA19SLA1,PCA19SLA1$med > 8)
PCA19SLA1LowBKMedOK %<>% mutate(med = 6.35)
PCA19SLA1LowBKMedOK %<>% mutate(mednorm = (log2-med))
PCA19SLA1LowBK = subset(PCA19SLA1,PCA19SLA1$med > 8)
PCA19SLA1 = anti_join(PCA19SLA1,PCA19SLA1LowBK)
PCA19SLA1 = rbind(PCA19SLA1,PCA19SLA1LowBKMedOK)

# Combination of the replicates.

PCA19SLA1MedOKRepCombined = PCA19SLA1 %>% group_by(DHFR1.2, DHFR3) %>%
  summarise(n(),medsize = median(size,na.rm=T), medrep = median(mednorm,na.rm=T))

# Distribution of all medrep, which corresponds to the PCA scores of the PPIs. Based on the distribution, most false PPIs have a medrep under 1.22 (above threshold = 3.5th percentile).

disMednormFinal = ggplot(data=PCA19SLA1MedOKRepCombined, aes(x = medrep, colour = DHFR1.2)) +geom_density() + geom_vline(xintercept=1.22)
pdf("distributionMedrep_1219Analysis.pdf",7,10)
disMednormFinal + theme(legend.position="none")
dev.off()

# However, there are many weak PPIs that are different between the four Sla1 WT baits (WT, WT*|WT|WT, WT|WT*|WT and WT|WT|WT*). These PPIs might be false PPIs. Increasing the threshold to 1.8 remove most of the differences in PPIs between the Sla1 WT baits.
# This threshold correspond to the top 2.75% of all detected colonies.

quantile(PCA19SLA1MedOKRepCombined$medrep,.9725)

disMednormFinal = ggplot(data=PCA19SLA1MedOKRepCombined, aes(x = medrep, colour = DHFR1.2)) +geom_density() + geom_vline(xintercept=1.80)
pdf("distributionMedrepFinal_1219Analysis.pdf",7,10)
disMednormFinal + theme(legend.position="none")
dev.off()

PCA19SLA1MedOKRepCombinedINTFInal= subset(PCA19SLA1MedOKRepCombined, PCA19SLA1MedOKRepCombined$medrep >= 1.80)
write.csv(PCA19SLA1MedOKRepCombinedINTFInal,"PCAInteractionsFinal_1219Analysis.csv")

# Remove all N = 1. We only consider PPIs with at least 2 replicates left after removing the colonies that didn't grow in the pipeline.

PCA19SLA1MedOKRepCombinedINTFInal %<>% filter(`n()`!=1)
write.csv(PCA19SLA1MedOKRepCombinedINTFInal,"PCAInteractionsN2Final_1219AnalysisN2.csv")
# This dataframe correspond to the analyzed data in the dataset that is available with the article. Some column names and bait names have been changed for the paper.

# Heatmap Figure S5A.

Baits =as.list(unique(PCA19SLA1MedOKRepCombinedINTFInal$DHFR1.2))
Preys = as.list(unique(PCA19SLA1MedOKRepCombinedINTFInal$DHFR3))
heatmapdf = data.frame(matrix(vector(), 33, 0, dimnames = list(c(Baits), c())))
for (i in Preys) {
  subdf = subset(PCA19SLA1MedOKRepCombinedINTFInal, PCA19SLA1MedOKRepCombinedINTFInal$DHFR3 == i)
  for (j in Baits) {
    if(j %in% subdf$DHFR1.2) { subsubdf = subset(subdf, subdf$DHFR1.2 == j)
    heatmapdf[j, i] = subsubdf$medrep }
    else { heatmapdf[j, i] = 0 }
  }
}
heatmap = data.matrix(heatmapdf, rownames.force = TRUE)
pdf("HeatmapSLA1Allinteractions_0120Analysis.pdf",12, 12)
heatmap.2(heatmap,dendrogram = "row",cexRow = 0.7, cexCol = 0.4,col=mycol, scale ="none", trace = "none", density.info = "none",margins = c(5,20))
dev.off()

# Heatmap Abp1 SH3 in Sla1, Figure 2E

SLA1ABP1Stuffed = subset(PCA19SLA1MedOKRepCombinedINTFInal, DHFR1.2 %in% c("SLA1 WT|WT|WT","SLA1 ABP1|WT|WT","SLA1 WT|ABP1|WT","SLA1 WT|WT|ABP1"))
Baits =as.list(unique(SLA1ABP1Stuffed$DHFR1.2))
Preys = as.list(unique(SLA1ABP1Stuffed$DHFR3))
heatmapdf = data.frame(matrix(vector(), 4, 0, dimnames = list(c(Baits), c())))
for (i in Preys) {
  subdf = subset(SLA1ABP1Stuffed, SLA1ABP1Stuffed$DHFR3 == i)
  for (j in Baits) {
    if(j %in% subdf$DHFR1.2) { subsubdf = subset(subdf, subdf$DHFR1.2 == j)
    heatmapdf[j, i] = subsubdf$medrep }
    else { heatmapdf[j, i] = 0 }
  }
}
heatmap = data.matrix(heatmapdf, rownames.force = TRUE)
pdf("HeatmapSLA1ABP1SWAPInteractions_0120Analysis.pdf",12, 10)
heatmap.2(heatmap,dendrogram = "row",cexRow = 2.0, cexCol = 0.4,col=mycol, scale ="none", trace = "none", density.info = "none",margins = c(5,20))
dev.off()

