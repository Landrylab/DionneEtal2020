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
library(GGally)

# Script for the analysis of the first PCA screen, including WT, SH3-deleted and Abp1 SH3-swapped baits.

# Gitter Analysis of plate pictures.

Ref="C:/Users/ugols/Desktop/PCA_Oct2018/Gitter_Analysis_S2_MTX2Last/d000097_300_001_18-11-10_08-39-46.JPG"
Analyse="C:/Users/ugols/Desktop/PCA_Oct2018/Gitter_Analysis_S2_MTX2Last"
path_dat="C:/Users/ugols/Desktop/PCA_Oct2018/Gitter_Analysis_S2_MTX2Last"
path_grid="C:/Users/ugols/Desktop/PCA_Oct2018/Gitter_Analysis_S2_MTX2Last"
gitter.batch(image.files=Analyse, ref.image.file=Ref, plate.format=1536, verbose='p',grid.save=path_grid, dat.save=path_dat)


# Merge Plate definitions, Prey positions on array with the Gitter analysis.

platesDHFR12_def_Oct2018 <- read.csv("plates_bait_UD_Oct2018_VF.csv", header = T, sep = ";")
platesDHFR3_def_Oct2018 <- read.csv("PreysArray_Oct2018.csv", header = T, sep = ",")

df = data.frame(matrix(vector(), 0, 10, dimnames = list(c(), c("row", "col", "size", "circularity", "flags", "plate", "DHFR1.2", "condition", "rep", "DHFR3"))))
head(df)
for (i in 1:nrow(platesDHFR12_def_Oct2018)) {
  data = read.table(paste("Gitter_Analysis_S2_MTX2Last/",as.character(platesDHFR12_def_Oct2018[i,"dat_Gitter"]), sep=""), sep="\t", header=F, stringsAsFactor=F, strip.white = T,comment.char = "#", na.strings="")
  names(data) <- c("row", "col", "size", "circularity", "flags")
  plate = platesDHFR12_def_Oct2018[i, "plate"]
  DHFR1.2 = platesDHFR12_def_Oct2018[i, "DHFR1.2"]
  condition = platesDHFR12_def_Oct2018[i, "condition"]
  rep = platesDHFR12_def_Oct2018[i, "rep"]
  subdf <- cbind(data, rep(plate, nrow(data)), rep(DHFR1.2, nrow(data)), rep(condition, nrow(data)), rep(rep, nrow(data)))
  subdf <- merge(subdf, platesDHFR3_def_Oct2018[,c("col", "row", "orf")], by.x=c("col", "row"), by.y=c("col", "row"))
  names(subdf) <- c("row", "col", "size", "circularity", "flags", "plate", "DHFR1.2", "condition", "rep", "DHFR3")
  df = rbind(df, subdf)
}
write.csv(df,"AlldataOct2018.csv")

# The datasets available with the article that are corresponding to the raw data for the first DHFR-PCA screen are from this dataset. Some colunm names and some baits names have been edited manually for the paper.

df$size = as.numeric(as.character(df$size))
df$col = as.numeric(as.character(df$col))
df$row = as.numeric(as.character(df$row))
df$plate = as.numeric(as.character(df$plate))
df$rep = as.numeric(as.character(df$rep))
df$DHFR3 = as.character(df$DHFR3)
df$DHFR1.2 = as.character(df$DHFR1.2)

# Analysis of the colonies from the diploid selection S2. Remove colonies that didn't grow after the diploid selection and then remove the S2 data.

S2= subset(df, df$condition == "S2")
S2null = subset(df, df$condition == "S2" & df$size < 600)
df = anti_join(df, S2null, by = c("row", "col", "plate"))
df %<>% filter(condition!="S2")

# Log2 transformation of colony size and remove positions with no growth in the DHFR-PCA selection.

df %<>% mutate(adjust = size+1)
df %<>% mutate(log2 = log2(adjust))
df %<>% mutate(baitXreplicate = interaction(DHFR1.2,rep))
df = subset(df, df$log2 > 2)

# Normalize the colony size to plate background.

df %<>% 
  group_by(baitXreplicate) %>%
  mutate(med = median(log2, na.rm=T)) %>%
  ungroup()
df %<>% mutate(mednorm = log2 - med)
write.csv(df,"AlldataBeforeCombining_Oct2018.csv")
AlldataNorm = read.csv("AlldataBeforeCombining_Oct2018.csv")
AlldataNorm$DHFR1.2 <- (as.character(AlldataNorm$DHFR1.2))
AlldataNorm$DHFR3 <- (as.character(AlldataNorm$DHFR3))

# Remove Border colonies.

AlldataNorm %<>% filter(DHFR3!="border")

# Look at the distribution of plate background.

DistributionMedian = ggplot(data = AlldataNorm, aes(x = med, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(0,15,0.5)) + geom_density() + geom_vline(xintercept=0)
pdf("DistributionMedian_1219Final.pdf")
DistributionMedian + theme(legend.position = "none")
dev.off()

# 2 plates have almost no background level which results in high plate median (colonies with no growth were removed).
# The 2 plates are Sla1-1 in Abp1, they are assigned the plate median corresponding to the average of Abp1 WT and Abp1 in Abp1 plates (7.13, 7.16, 7.2, 6.3). The assigned plate median is 7.0.

SLA11inABP1OK = subset(AlldataNorm, AlldataNorm$DHFR1.2 =="SLA1-1 in ABP1")
SLA11inABP1OK %<>% mutate(med = 7.0)
SLA11inABP1OK %<>% mutate(mednorm = (log2 - med))
SLA11inABP1 = subset(AlldataNorm, AlldataNorm$DHFR1.2 =="SLA1-1 in ABP1")

AlldataNorm = anti_join(AlldataNorm,SLA11inABP1)
AlldataNorm = rbind(AlldataNorm,SLA11inABP1OK)

# Combination of the replicates.

AlldataNormRepComb = AlldataNorm %>% group_by(DHFR1.2, DHFR3) %>%
  summarise(n(),medsize = median(size,na.rm=T), medrep = median(mednorm,na.rm=T))

# Distribution of all medrep, which corresponds to the PCA scores of the PPIs. Based on the distribution, most false PPIs have a medrep under 1.03 (above the threshold = 8th percentile).

DistributionRepComb2 = ggplot(data = AlldataNormRepComb, aes(x = medrep, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(-10,7,0.5)) + geom_density() + geom_vline(xintercept=1.03)
pdf("DistributionRepComb2_1219Final.pdf")
DistributionRepComb2 + theme(legend.position = "none")
dev.off()

# However, there are many weak PPIs that are different between Abp1 WT and the Abp1 in Abp1 control. These PPIs might be false PPIs. Increasing the threshold to 1.3 remove most of the differences in PPIs between Abp1 WT and the control.
# 1.3 medrep (PCA score) corresponds to the top 7.3% of all detected colonies.

DistributionRepCombFinal = ggplot(data = AlldataNormRepComb, aes(x = medrep, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(-10,7,0.5)) + geom_density() + geom_vline(xintercept=1.3)
pdf("DistributionRepCombFinal_1219Final.pdf")
DistributionRepCombFinal + theme(legend.position = "none")
dev.off()
AlldataNormRepCombINT= subset(AlldataNormRepComb, AlldataNormRepComb$medrep >= 1.3)

#Remove all N = 1 (very few PPIs, around 20 total). We only consider PPIs with at least 2 replicates left after removing the colonies that didn't grow in the pipeline.

AlldataNormRepCombINT %<>% filter(`n()`!=1) 



## Analysis of PPIs that are affected by SH3 deletions by comparing the PPIs of WT and SH3-deleted (stuffed) baits.

SH3ProtsInt = subset(AlldataNormRepCombINT, DHFR1.2 %in% c("ABP1","BBC1","BEM1","BOI1","BOI2","BUD14","BZZ1","CDC25","CYK3","FUS1","HOF1","HSE1","LSB1","LSB3","MYO3","MYO5","NBP2","PEX13","PIN3","RVS167","SHO1","SLA1"))
SH3StuffedInt = subset(AlldataNormRepCombINT, DHFR1.2 %in% c("ABP1 Stuffed","BBC1 Stuffed","BEM1-1 Stuffed","BEM1-2 Stuffed","BEM1 doubled Stuffed","BOI1 Stuffed","BOI2 Stuffed","BUD14 Stuffed","BZZ1-1 Stuffed","BZZ1-2 Stuffed","BZZ1 doubled Stuffed","CDC25 Stuffed","CYK3 Stuffed","FUS1 Stuffed","HOF1 Stuffed","HSE1 Stuffed","LSB1 Stuffed","LSB3 Stuffed","MYO3 Stuffed","MYO5 Stuffed","NBP2 Stuffed","PEX13 Stuffed","PIN3 Stuffed","RVS167 Stuffed","SHO1 Stuffed","SLA1-1 Stuffed","SLA1-2 Stuffed","SLA1-3 Stuffed","SLA1 triple Stuffed"))

write.csv(SH3ProtsInt,"InteractionsPCAOct2018SH3WT_1219Final.csv")
write.csv(SH3StuffedInt,"InteractionsPCAOct2018Suffed_1219Final.csv")

# Manually removed the stuffed from the DHFR1.2 protein names to compare the dataframes.
SH3StuffedInt = read.csv("InteractionsPCAOct2018Suffed_1219FinalForanalysis.csv")
SH3StuffedInt$DHFR1.2 <- (as.character(SH3StuffedInt$DHFR1.2))
SH3StuffedInt$DHFR3 <- (as.character(SH3StuffedInt$DHFR3))

# Interactions that are lost or gained upon SH3 domain removal.
# Prepare the dataframes to compare the PPIs of WT and SH3-deleted (stuffed) baits.

BEM1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BEM1")
BEM1 %<>% ungroup()
BEM1%<>% mutate(DHFR1.2 = "BEM1-1")
BEM1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,BEM1)

BEM1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BEM1")
BEM1 %<>% ungroup()
BEM1%<>% mutate(DHFR1.2 = "BEM1-2")
BEM1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,BEM1)

BEM1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BEM1")
SH3ProtsInt = anti_join(SH3ProtsInt,BEM1)

BZZ1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BZZ1")
BZZ1 %<>% ungroup()
BZZ1%<>% mutate(DHFR1.2 = "BZZ1-1")
BZZ1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,BZZ1)

BZZ1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BZZ1")
BZZ1 %<>% ungroup()
BZZ1%<>% mutate(DHFR1.2 = "BZZ1-2")
BZZ1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,BZZ1)

BZZ1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "BZZ1")
SH3ProtsInt = anti_join(SH3ProtsInt,BZZ1)

SLA1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "SLA1")
SLA1 %<>% ungroup()
SLA1%<>% mutate(DHFR1.2 = "SLA1-1")
SLA1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,SLA1)

SLA1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "SLA1")
SLA1 %<>% ungroup()
SLA1%<>% mutate(DHFR1.2 = "SLA1-2")
SLA1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,SLA1)

SLA1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "SLA1")
SLA1 %<>% ungroup()
SLA1%<>% mutate(DHFR1.2 = "SLA1-3")
SLA1%<>% group_by(DHFR1.2)
SH3ProtsInt = rbind(SH3ProtsInt,SLA1)

SLA1 = subset(SH3ProtsInt,SH3ProtsInt$DHFR1.2 == "SLA1")
SH3ProtsInt = anti_join(SH3ProtsInt,SLA1)

SH3depIntLost = anti_join(SH3ProtsInt,SH3StuffedInt, by = c("DHFR1.2","DHFR3"))

SH3depIntGained = anti_join(SH3StuffedInt,SH3ProtsInt, by = c("DHFR1.2","DHFR3"))
SH3depIntGained%<>% filter(DHFR1.2!="BEM1 doubled")
SH3depIntGained%<>% filter(DHFR1.2!="BZZ1 doubled")
SH3depIntGained%<>% filter(DHFR1.2!="SLA1 triple")

# Comparison of the PCA scores (medrep) of the SH3-deleted (stuffed) baits with the WT baits (calculation of a ratio of PCA scores).

IntNotAffected = anti_join(SH3ProtsInt,SH3depIntLost, by = c("DHFR1.2","DHFR3"))
IntNotAffected%<>% mutate(Condition = "WT")
SH3StuffedInt%<>% mutate(Condition = "STUFFED")

RatioCalulationStuffer = merge(IntNotAffected,SH3StuffedInt, by = c("DHFR1.2","DHFR3"))
RatioCalulationStuffer%<>%mutate(Ratio = (medrep.y/medrep.x))

# Distribution of the Ratios.

DistributionRatioStuffer = ggplot(data = RatioCalulationStuffer, aes(x = Ratio, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(-10,7,0.5)) + geom_density() + geom_vline(xintercept=0)
pdf("DistributionRatioStuffer_1219Final.pdf")
DistributionRatioStuffer + theme(legend.position = "none")
dev.off()

# The cutoff for negatively or positively affected interactions is the same percentile as the cutoff for the interactions (bottom 7.3% = Ratio 0.58 and top 7.3% = Ratio 1.31).

quantile(RatioCalulationStuffer$Ratio, .927)
quantile(RatioCalulationStuffer$Ratio, .073)

# Distribution of the Ratios with the thresholds.

DistributionRatioStuffer = ggplot(data = RatioCalulationStuffer, aes(x = Ratio, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(-10,7,0.5)) + geom_density() + geom_vline(xintercept=c(0,0.58,1.31))
pdf("DistributionRatioStuffer_1219Final.pdf")
DistributionRatioStuffer + theme(legend.position = "none")
dev.off()

PositivelyAffectedInt = subset(RatioCalulationStuffer,RatioCalulationStuffer$Ratio >=1.31)
PositivelyAffectedInt[,c(3:6)]=NULL
PositivelyAffectedInt[,6]=NULL
colnames(PositivelyAffectedInt) = c("DHFR1.2","DHFR3","n()","medsize","medrep","Ratio")

NegativelyAffectedInt = subset(RatioCalulationStuffer,RatioCalulationStuffer$Ratio <= 0.58)
NegativelyAffectedInt[,c(3:6)]=NULL
NegativelyAffectedInt[,6]=NULL
colnames(NegativelyAffectedInt) = c("DHFR1.2","DHFR3","n()","medsize","medrep","Ratio")

# Combine everything in one dataframe. Ratio of 0.1 for the lost interactions and 2.5 for the gained. Score: 0 = SH3-independent PPIs, 1 = Lost or negatively affected, 2 = Gained or positively affected.

SH3depIntLost%<>% mutate(Ratio = 0.1)
SH3depIntLost%<>% mutate(Score = 1)
SH3depIntLost%<>% group_by(DHFR1.2)
SH3depIntGained%<>% mutate(Ratio = 2.5)
SH3depIntGained%<>% mutate(Score = 2)
SH3depIntGained%<>% group_by(DHFR1.2)
colnames(SH3depIntGained) = c("DHFR1.2","DHFR3","n()","medsize","medrep","Ratio","Score")

SH3depInteractions = rbind(SH3depIntLost,SH3depIntGained)

PositivelyAffectedInt%<>% mutate(Score = 2)
PositivelyAffectedInt%<>% group_by(DHFR1.2)
SH3depInteractions = rbind(SH3depInteractions,PositivelyAffectedInt)

NegativelyAffectedInt%<>% mutate(Score = 1)
NegativelyAffectedInt%<>% group_by(DHFR1.2)
SH3depInteractions = rbind(SH3depInteractions,NegativelyAffectedInt)

SH3IntNegAffec = rbind(SH3depIntLost,NegativelyAffectedInt)
SH3IntPosAffec = rbind(SH3depIntGained,PositivelyAffectedInt)

SH3IndepInteractions = anti_join(SH3ProtsInt,SH3depInteractions, by = c("DHFR1.2","DHFR3"))

SH3IndepInteractions%<>% mutate(Ratio = NA)
SH3IndepInteractions%<>% mutate(Score = 0)
SH3IndepInteractions%<>%group_by(DHFR1.2)

AllSH3Interactions = rbind(SH3depInteractions,SH3IndepInteractions)

write.csv(AllSH3Interactions,"PCAOct18_AllSH3interactions_0120Analysis.csv")
# This corresponds to the analyzed dataset available with the article. Some colunm names and some bait names have been edited manually for the paper.

## Heatmap generation for Figure 1B.

AllSH3Interactions = read.csv("PCAOct18_AllSH3interactions_0120Analysis.csv", header =T)
AllSH3Interactions$DHFR1.2 <- (as.character(AllSH3Interactions$DHFR1.2))
AllSH3Interactions$DHFR3 <- (as.character(AllSH3Interactions$DHFR3))

#Heatmap. SH3-independent PPIs (grey), SH3-dependent (blue), SH3-inhibited (red) and no PPIs (white). Colors from RdBu.

my_palette = colorRampPalette(c("#DCDCDC", "#2166AC", "#B2182B","#F8F8F8"))(n = 4)

Baits =as.list(unique(AllSH3Interactions$DHFR1.2))
Preys = as.list(unique(AllSH3Interactions$DHFR3))
heatmapdf = data.frame(matrix(vector(), 26, 0, dimnames = list(c(Baits), c())))
for (i in Preys) {
  subdf = subset(AllSH3Interactions, AllSH3Interactions$DHFR3 == i)
  for (j in Baits) {
    if(j %in% subdf$DHFR1.2) { subsubdf = subset(subdf, subdf$DHFR1.2 == j)
    heatmapdf[j, i] = subsubdf$Score }
    else { heatmapdf[j, i] = 3 }
  }
}
heatmap = data.matrix(heatmapdf, rownames.force = TRUE)
pdf("HeatmapPCAOCTAllSH3INT_0120Analysis_Fig0320.pdf",7, 10)
heatmap.2(heatmap,dendrogram = "row",cexRow = 0.75, cexCol = 0.25,col=my_palette, scale ="none", trace = "none", density.info = "none",margins = c(5,7))
dev.off()
# The final figure has also been edited in illustrator.

#### Determine the correlation between the different replicates of PPIs above the treshold for the WT and SH3 stuffed baits.

# Get PPIs that are above the treshold.

SH3WTandStuffedPPIs = subset(AlldataNormRepCombINT, DHFR1.2 %in% c("ABP1","BBC1","BEM1","BOI1","BOI2","BUD14","BZZ1","CDC25","CYK3","FUS1","HOF1","HSE1","LSB1","LSB3","MYO3","MYO5","NBP2","PEX13","PIN3","RVS167","SHO1","SLA1","ABP1 Stuffed","BBC1 Stuffed","BEM1-1 Stuffed","BEM1-2 Stuffed","BEM1 doubled Stuffed","BOI1 Stuffed","BOI2 Stuffed","BUD14 Stuffed","BZZ1-1 Stuffed","BZZ1-2 Stuffed","BZZ1 doubled Stuffed","CDC25 Stuffed","CYK3 Stuffed","FUS1 Stuffed","HOF1 Stuffed","HSE1 Stuffed","LSB1 Stuffed","LSB3 Stuffed","MYO3 Stuffed","MYO5 Stuffed","NBP2 Stuffed","PEX13 Stuffed","PIN3 Stuffed","RVS167 Stuffed","SHO1 Stuffed","SLA1-1 Stuffed","SLA1-2 Stuffed","SLA1-3 Stuffed","SLA1 triple Stuffed"))
SH3WTandStuffedPPIs %<>% mutate(BaitxPrey = paste(DHFR1.2, DHFR3, sep="_"))

# Get PPIs replicates with normalized values.

SH3WTandStuffedPPIsReplicates = subset(AlldataNorm, DHFR1.2 %in% c("ABP1","BBC1","BEM1","BOI1","BOI2","BUD14","BZZ1","CDC25","CYK3","FUS1","HOF1","HSE1","LSB1","LSB3","MYO3","MYO5","NBP2","PEX13","PIN3","RVS167","SHO1","SLA1","ABP1 Stuffed","BBC1 Stuffed","BEM1-1 Stuffed","BEM1-2 Stuffed","BEM1 doubled Stuffed","BOI1 Stuffed","BOI2 Stuffed","BUD14 Stuffed","BZZ1-1 Stuffed","BZZ1-2 Stuffed","BZZ1 doubled Stuffed","CDC25 Stuffed","CYK3 Stuffed","FUS1 Stuffed","HOF1 Stuffed","HSE1 Stuffed","LSB1 Stuffed","LSB3 Stuffed","MYO3 Stuffed","MYO5 Stuffed","NBP2 Stuffed","PEX13 Stuffed","PIN3 Stuffed","RVS167 Stuffed","SHO1 Stuffed","SLA1-1 Stuffed","SLA1-2 Stuffed","SLA1-3 Stuffed","SLA1 triple Stuffed"))
SH3WTandStuffedPPIsReplicates %<>% mutate(BaitxPrey = paste(DHFR1.2, DHFR3, sep="_"))

#Merge both.

ReplicatesCorrelation = merge(SH3WTandStuffedPPIs,SH3WTandStuffedPPIsReplicates, by = "BaitxPrey") 
ReplicatesCorrelationDF = subset(ReplicatesCorrelation, select=c(BaitxPrey,size,medsize,log2,mednorm,medrep))

ReplicatesCorrelationDFFInal = subset(ReplicatesCorrelationDF, select = c(BaitxPrey,mednorm))

agg = aggregate (ReplicatesCorrelationDFFInal,by = list(ReplicatesCorrelationDFFInal$BaitxPrey), FUN = unique)

# Separate the values in separate columns.

setDT(agg)[, paste0("mednorm", 1:8) := tstrsplit(mednorm, ",")]
setDT(agg)[, paste0("mednorm1", 1:2) := tstrsplit(mednorm1, "(", type.convert = TRUE, fixed = TRUE)]
setDT(agg)[, paste0("mednorm4", 1:2) := tstrsplit(mednorm4, ")", type.convert = TRUE, fixed = TRUE)]
setDT(agg)[, paste0("mednorm2", 1:2) := tstrsplit(mednorm2, ")", type.convert = TRUE, fixed = TRUE)]
setDT(agg)[, paste0("mednorm3", 1:2) := tstrsplit(mednorm3, ")", type.convert = TRUE, fixed = TRUE)]

ComparisonReplicatesPPIs = subset(agg, select=c(BaitxPrey,mednorm12,mednorm21,mednorm31,mednorm41))
colnames(ComparisonReplicatesPPIs) = c("BaitxPrey","Rep1","Rep2","Rep3","Rep4")

ComparisonReplicatesPPIs$Rep2 = as.numeric(as.character(ComparisonReplicatesPPIs$Rep2))
ComparisonReplicatesPPIs$Rep3 = as.numeric(as.character(ComparisonReplicatesPPIs$Rep3))

# correlation calculation and Heatmap generation.

pdf("CorrelationPPIsRep_0420.pdf",10, 10)
ggcorr(ComparisonReplicatesPPIs[, 2:5], low = "steelblue", mid = "white", high = "darkred")
CorrelationData = ggcorr(ComparisonReplicatesPPIs[, 2:5])$data
dev.off()


##### Abp1 SH3-swapped analysis

# This is the dataset that correspond to the analyzed data available with the paper. Datasets have been edited manually to change some column names and bait names for the paper.

ABP1Swap=AlldataNormRepCombINT[grep("ABP1", AlldataNormRepCombINT$DHFR1.2),]
write.csv(ABP1Swap,"PCAOct18ABP1SWAPInt_0120Analysis.csv")

# Heatmap generation for Figure 2B.
# Remove the in Abp1 nomenclature.

Abp1Swap = read.csv("PCAOct18ABP1SWAPInt_0120Analysis.csv")
Abp1Swap$DHFR1.2 <- (as.character(Abp1Swap$DHFR1.2))
Abp1Swap$DHFR3 <- (as.character(Abp1Swap$DHFR3))

Abp1SwapWT = subset(Abp1Swap,Abp1Swap$DHFR1.2 == "ABP1")
Abp1Swap = anti_join(Abp1Swap,Abp1SwapWT)
Abp1SwapWT %<>% ungroup()
Abp1SwapWT %<>% mutate(DHFR1.2 = "ABP1WT")
Abp1Swap = rbind(Abp1Swap, Abp1SwapWT)

Abp1SwapStuffed = subset(Abp1Swap,Abp1Swap$DHFR1.2 == "ABP1 Stuffed")
Abp1Swap = anti_join(Abp1Swap,Abp1SwapStuffed)
Abp1SwapStuffed %<>% ungroup()
Abp1SwapStuffed %<>% mutate(DHFR1.2 = "ABP1Stuffed")
Abp1Swap = rbind(Abp1Swap, Abp1SwapStuffed)

setDT(Abp1Swap)[, paste0("DHFR1.2", 1:3) := tstrsplit(DHFR1.2, " ")]
Abp1Swap = subset(Abp1Swap, select=c(DHFR1.21,2:5))
colnames(Abp1Swap) = c("DHFR1.2","DHFR3","n","medsize","medrep")

#Add number of interactions per prey to order the heatmap.

Abp1Swap = add_count(Abp1Swap, DHFR3)
Abp1Swap = Abp1Swap[order(-Abp1Swap$n,Abp1Swap$DHFR3),]

#Heatmap ordered by number of PPIs by preys.

Preys = as.list(unique(Abp1Swap$DHFR3))
Strains = as.list(unique(Abp1Swap$DHFR1.2))
heatmapdf = data.frame(matrix(vector(), 33, 0, dimnames = list(c(Strains), c())))
for (i in Preys) {
  subdf = subset(Abp1Swap, Abp1Swap$DHFR3 == i)
  for (j in Strains) {
    if(j %in% subdf$DHFR1.2) { subsubdf = subset(subdf, subdf$DHFR1.2 == j)
    heatmapdf[j, i] = subsubdf$medrep }
    else { heatmapdf[j, i] = 0 }
  }
}
heatmap = data.matrix(heatmapdf, rownames.force = TRUE)
pdf("HeatmapPCAABP1SWAP_0120Analysis_Fig0320.pdf",7, 10)
heatmap.2(heatmap,dendrogram = "row",Colv = FALSE,cexRow = 0.6, cexCol = 0.35,col=mycol, scale ="none", trace = "none", density.info = "none",margins = c(5,7))
dev.off()
# The final figure has also been edited in illustrator.

#### Determination of PPIs that are lost or gained following domain swapping in Abp1 when compared to Abp1 in Abp1.

ABP1Swap=AlldataNormRepCombINT[grep("ABP1", AlldataNormRepCombINT$DHFR1.2),]

ABP1inABP1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "ABP1 in ABP1")
BBC1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BBC1 in ABP1")
BBC1Lost = anti_join(ABP1inABP1,BBC1, by = "DHFR3")
BBC1Lost %<>% ungroup()
BBC1Lost%<>% mutate(DHFR1.2 = "BBC1 in ABP1")
BBC1Lost%<>% group_by(DHFR1.2)
BBC1Ratio = anti_join(ABP1inABP1,BBC1Lost, by = "DHFR3")
BBC1Ratio = merge(BBC1,BBC1Ratio, by = "DHFR3")
BBC1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BBC1Ratio[,c(6:9)]=NULL
colnames(BBC1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BBC1gain = anti_join(BBC1, ABP1inABP1, by ="DHFR3")

BEM11 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BEM1-1 in ABP1")
BEM11Lost = anti_join(ABP1inABP1,BEM11, by = "DHFR3")
BEM11Lost %<>% ungroup()
BEM11Lost%<>% mutate(DHFR1.2 = "BEM1-1 in ABP1")
BEM11Lost%<>% group_by(DHFR1.2)
BEM11Ratio = anti_join(ABP1inABP1,BEM11Lost, by = "DHFR3")
BEM11Ratio = merge(BEM11,BEM11Ratio, by = "DHFR3")
BEM11Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BEM11Ratio[,c(6:9)]=NULL
colnames(BEM11Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BEM11gain = anti_join(BEM11, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(BBC1Lost,BEM11Lost)
ABP1SWAPRatio = rbind(BBC1Ratio, BEM11Ratio)
ABP1SWAPGain = rbind(BBC1gain,BEM11gain)

BEM12 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BEM1-2 in ABP1")
BEM12Lost = anti_join(ABP1inABP1,BEM12, by = "DHFR3")
BEM12Lost %<>% ungroup()
BEM12Lost%<>% mutate(DHFR1.2 = "BEM1-2 in ABP1")
BEM12Lost%<>% group_by(DHFR1.2)
BEM12Ratio = anti_join(ABP1inABP1,BEM12Lost, by = "DHFR3")
BEM12Ratio = merge(BEM12,BEM12Ratio, by = "DHFR3")
BEM12Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BEM12Ratio[,c(6:9)]=NULL
colnames(BEM12Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BEM12gain = anti_join(BEM12, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BEM12Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BEM12Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BEM12gain)

BOI1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BOI1 in ABP1")
BOI1Lost = anti_join(ABP1inABP1,BOI1, by = "DHFR3")
BOI1Lost %<>% ungroup()
BOI1Lost%<>% mutate(DHFR1.2 = "BOI1 in ABP1")
BOI1Lost%<>% group_by(DHFR1.2)
BOI1Ratio = anti_join(ABP1inABP1,BOI1Lost, by = "DHFR3")
BOI1Ratio = merge(BOI1,BOI1Ratio, by = "DHFR3")
BOI1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BOI1Ratio[,c(6:9)]=NULL
colnames(BOI1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BOI1gain = anti_join(BOI1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BOI1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BOI1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BOI1gain)

BOI2 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BOI2 in ABP1")
BOI2Lost = anti_join(ABP1inABP1,BOI2, by = "DHFR3")
BOI2Lost %<>% ungroup()
BOI2Lost%<>% mutate(DHFR1.2 = "BOI2 in ABP1")
BOI2Lost%<>% group_by(DHFR1.2)
BOI2Ratio = anti_join(ABP1inABP1,BOI2Lost, by = "DHFR3")
BOI2Ratio = merge(BOI2,BOI2Ratio, by = "DHFR3")
BOI2Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BOI2Ratio[,c(6:9)]=NULL
colnames(BOI2Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BOI2gain = anti_join(BOI2, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BOI2Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BOI2Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BOI2gain)

BUD14 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BUD14 in ABP1")
BUD14Lost = anti_join(ABP1inABP1,BUD14, by = "DHFR3")
BUD14Lost %<>% ungroup()
BUD14Lost%<>% mutate(DHFR1.2 = "BUD14 in ABP1")
BUD14Lost%<>% group_by(DHFR1.2)
BUD14Ratio = anti_join(ABP1inABP1,BUD14Lost, by = "DHFR3")
BUD14Ratio = merge(BUD14,BUD14Ratio, by = "DHFR3")
BUD14Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BUD14Ratio[,c(6:9)]=NULL
colnames(BUD14Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BUD14gain = anti_join(BUD14, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BUD14Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BUD14Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BUD14gain)

BZZ11 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BZZ1-1 in ABP1")
BZZ11Lost = anti_join(ABP1inABP1,BZZ11, by = "DHFR3")
BZZ11Lost %<>% ungroup()
BZZ11Lost%<>% mutate(DHFR1.2 = "BZZ1-1 in ABP1")
BZZ11Lost%<>% group_by(DHFR1.2)
BZZ11Ratio = anti_join(ABP1inABP1,BZZ11Lost, by = "DHFR3")
BZZ11Ratio = merge(BZZ11,BZZ11Ratio, by = "DHFR3")
BZZ11Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BZZ11Ratio[,c(6:9)]=NULL
colnames(BZZ11Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BZZ11gain = anti_join(BZZ11, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BZZ11Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BZZ11Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BZZ11gain)

BZZ12 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "BZZ1-2 in ABP1")
BZZ12Lost = anti_join(ABP1inABP1,BZZ12, by = "DHFR3")
BZZ12Lost %<>% ungroup()
BZZ12Lost%<>% mutate(DHFR1.2 = "BZZ1-2 in ABP1")
BZZ12Lost%<>% group_by(DHFR1.2)
BZZ12Ratio = anti_join(ABP1inABP1,BZZ12Lost, by = "DHFR3")
BZZ12Ratio = merge(BZZ12,BZZ12Ratio, by = "DHFR3")
BZZ12Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
BZZ12Ratio[,c(6:9)]=NULL
colnames(BZZ12Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
BZZ12gain = anti_join(BZZ12, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,BZZ12Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, BZZ12Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,BZZ12gain)

CDC25 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "CDC25 in ABP1")
CDC25Lost = anti_join(ABP1inABP1,CDC25, by = "DHFR3")
CDC25Lost %<>% ungroup()
CDC25Lost%<>% mutate(DHFR1.2 = "CDC25 in ABP1")
CDC25Lost%<>% group_by(DHFR1.2)
CDC25Ratio = anti_join(ABP1inABP1,CDC25Lost, by = "DHFR3")
CDC25Ratio = merge(CDC25,CDC25Ratio, by = "DHFR3")
CDC25Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
CDC25Ratio[,c(6:9)]=NULL
colnames(CDC25Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
CDC25gain = anti_join(CDC25, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,CDC25Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, CDC25Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,CDC25gain)

CTTN = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "CTTN in ABP1")
CTTNLost = anti_join(ABP1inABP1,CTTN, by = "DHFR3")
CTTNLost %<>% ungroup()
CTTNLost%<>% mutate(DHFR1.2 = "CTTN in ABP1")
CTTNLost%<>% group_by(DHFR1.2)
CTTNRatio = anti_join(ABP1inABP1,CTTNLost, by = "DHFR3")
CTTNRatio = merge(CTTN,CTTNRatio, by = "DHFR3")
CTTNRatio %<>% mutate(Ratio = medrep.x/medrep.y)
CTTNRatio[,c(6:9)]=NULL
colnames(CTTNRatio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
CTTNgain = anti_join(CTTN, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,CTTNLost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, CTTNRatio)
ABP1SWAPGain = rbind(ABP1SWAPGain,CTTNgain)

CYK3 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "CYK3 in ABP1")
CYK3Lost = anti_join(ABP1inABP1,CYK3, by = "DHFR3")
CYK3Lost %<>% ungroup()
CYK3Lost%<>% mutate(DHFR1.2 = "CYK3 in ABP1")
CYK3Lost%<>% group_by(DHFR1.2)
CYK3Ratio = anti_join(ABP1inABP1,CYK3Lost, by = "DHFR3")
CYK3Ratio = merge(CYK3,CYK3Ratio, by = "DHFR3")
CYK3Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
CYK3Ratio[,c(6:9)]=NULL
colnames(CYK3Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
CYK3gain = anti_join(CYK3, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,CYK3Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, CYK3Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,CYK3gain)

DBNL = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "DBNL in ABP1")
DBNLLost = anti_join(ABP1inABP1,DBNL, by = "DHFR3")
DBNLLost %<>% ungroup()
DBNLLost%<>% mutate(DHFR1.2 = "DBNL in ABP1")
DBNLLost%<>% group_by(DHFR1.2)
DBNLRatio = anti_join(ABP1inABP1,DBNLLost, by = "DHFR3")
DBNLRatio = merge(DBNL,DBNLRatio, by = "DHFR3")
DBNLRatio %<>% mutate(Ratio = medrep.x/medrep.y)
DBNLRatio[,c(6:9)]=NULL
colnames(DBNLRatio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
DBNLgain = anti_join(DBNL, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,DBNLLost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, DBNLRatio)
ABP1SWAPGain = rbind(ABP1SWAPGain,DBNLgain)

FUS1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "FUS1 in ABP1")
FUS1Lost = anti_join(ABP1inABP1,FUS1, by = "DHFR3")
FUS1Lost %<>% ungroup()
FUS1Lost%<>% mutate(DHFR1.2 = "FUS1 in ABP1")
FUS1Lost%<>% group_by(DHFR1.2)
FUS1Ratio = anti_join(ABP1inABP1,FUS1Lost, by = "DHFR3")
FUS1Ratio = merge(FUS1,FUS1Ratio, by = "DHFR3")
FUS1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
FUS1Ratio[,c(6:9)]=NULL
colnames(FUS1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
FUS1gain = anti_join(FUS1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,FUS1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, FUS1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,FUS1gain)

HCLS1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "HCLS1 in ABP1")
HCLS1Lost = anti_join(ABP1inABP1,HCLS1, by = "DHFR3")
HCLS1Lost %<>% ungroup()
HCLS1Lost%<>% mutate(DHFR1.2 = "HCLS1 in ABP1")
HCLS1Lost%<>% group_by(DHFR1.2)
HCLS1Ratio = anti_join(ABP1inABP1,HCLS1Lost, by = "DHFR3")
HCLS1Ratio = merge(HCLS1,HCLS1Ratio, by = "DHFR3")
HCLS1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
HCLS1Ratio[,c(6:9)]=NULL
colnames(HCLS1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
HCLS1gain = anti_join(HCLS1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,HCLS1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, HCLS1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,HCLS1gain)

HOF1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "HOF1 in ABP1")
HOF1Lost = anti_join(ABP1inABP1,HOF1, by = "DHFR3")
HOF1Lost %<>% ungroup()
HOF1Lost%<>% mutate(DHFR1.2 = "HOF1 in ABP1")
HOF1Lost%<>% group_by(DHFR1.2)
HOF1Ratio = anti_join(ABP1inABP1,HOF1Lost, by = "DHFR3")
HOF1Ratio = merge(HOF1,HOF1Ratio, by = "DHFR3")
HOF1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
HOF1Ratio[,c(6:9)]=NULL
colnames(HOF1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
HOF1gain = anti_join(HOF1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,HOF1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, HOF1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,HOF1gain)

HSE1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "HSE1 in ABP1")
HSE1Lost = anti_join(ABP1inABP1,HSE1, by = "DHFR3")
HSE1Lost %<>% ungroup()
HSE1Lost%<>% mutate(DHFR1.2 = "HSE1 in ABP1")
HSE1Lost%<>% group_by(DHFR1.2)
HSE1Ratio = anti_join(ABP1inABP1,HSE1Lost, by = "DHFR3")
HSE1Ratio = merge(HSE1,HSE1Ratio, by = "DHFR3")
HSE1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
HSE1Ratio[,c(6:9)]=NULL
colnames(HSE1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
HSE1gain = anti_join(HSE1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,HSE1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, HSE1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,HSE1gain)

LSB1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "LSB1 in ABP1")
LSB1Lost = anti_join(ABP1inABP1,LSB1, by = "DHFR3")
LSB1Lost %<>% ungroup()
LSB1Lost%<>% mutate(DHFR1.2 = "LSB1 in ABP1")
LSB1Lost%<>% group_by(DHFR1.2)
LSB1Ratio = anti_join(ABP1inABP1,LSB1Lost, by = "DHFR3")
LSB1Ratio = merge(LSB1,LSB1Ratio, by = "DHFR3")
LSB1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
LSB1Ratio[,c(6:9)]=NULL
colnames(LSB1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
LSB1gain = anti_join(LSB1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,LSB1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, LSB1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,LSB1gain)

LSB3 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "LSB3 in ABP1")
LSB3Lost = anti_join(ABP1inABP1,LSB3, by = "DHFR3")
LSB3Lost %<>% ungroup()
LSB3Lost%<>% mutate(DHFR1.2 = "LSB3 in ABP1")
LSB3Lost%<>% group_by(DHFR1.2)
LSB3Ratio = anti_join(ABP1inABP1,LSB3Lost, by = "DHFR3")
LSB3Ratio = merge(LSB3,LSB3Ratio, by = "DHFR3")
LSB3Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
LSB3Ratio[,c(6:9)]=NULL
colnames(LSB3Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
LSB3gain = anti_join(LSB3, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,LSB3Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, LSB3Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,LSB3gain)

MYO3 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "MYO3 in ABP1")
MYO3Lost = anti_join(ABP1inABP1,MYO3, by = "DHFR3")
MYO3Lost %<>% ungroup()
MYO3Lost%<>% mutate(DHFR1.2 = "MYO3 in ABP1")
MYO3Lost%<>% group_by(DHFR1.2)
MYO3Ratio = anti_join(ABP1inABP1,MYO3Lost, by = "DHFR3")
MYO3Ratio = merge(MYO3,MYO3Ratio, by = "DHFR3")
MYO3Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
MYO3Ratio[,c(6:9)]=NULL
colnames(MYO3Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
MYO3gain = anti_join(MYO3, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,MYO3Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, MYO3Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,MYO3gain)

MYO5 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "MYO5 in ABP1")
MYO5Lost = anti_join(ABP1inABP1,MYO5, by = "DHFR3")
MYO5Lost %<>% ungroup()
MYO5Lost%<>% mutate(DHFR1.2 = "MYO5 in ABP1")
MYO5Lost%<>% group_by(DHFR1.2)
MYO5Ratio = anti_join(ABP1inABP1,MYO5Lost, by = "DHFR3")
MYO5Ratio = merge(MYO5,MYO5Ratio, by = "DHFR3")
MYO5Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
MYO5Ratio[,c(6:9)]=NULL
colnames(MYO5Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
MYO5gain = anti_join(MYO5, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,MYO5Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, MYO5Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,MYO5gain)

NBP2 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "NBP2 in ABP1")
NBP2Lost = anti_join(ABP1inABP1,NBP2, by = "DHFR3")
NBP2Lost %<>% ungroup()
NBP2Lost%<>% mutate(DHFR1.2 = "NBP2 in ABP1")
NBP2Lost%<>% group_by(DHFR1.2)
NBP2Ratio = anti_join(ABP1inABP1,NBP2Lost, by = "DHFR3")
NBP2Ratio = merge(NBP2,NBP2Ratio, by = "DHFR3")
NBP2Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
NBP2Ratio[,c(6:9)]=NULL
colnames(NBP2Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
NBP2gain = anti_join(NBP2, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,NBP2Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, NBP2Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,NBP2gain)

PEX13 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "PEX13 in ABP1")
PEX13Lost = anti_join(ABP1inABP1,PEX13, by = "DHFR3")
PEX13Lost %<>% ungroup()
PEX13Lost%<>% mutate(DHFR1.2 = "PEX13 in ABP1")
PEX13Lost%<>% group_by(DHFR1.2)
PEX13Ratio = anti_join(ABP1inABP1,PEX13Lost, by = "DHFR3")
PEX13Ratio = merge(PEX13,PEX13Ratio, by = "DHFR3")
PEX13Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
PEX13Ratio[,c(6:9)]=NULL
colnames(PEX13Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
PEX13gain = anti_join(PEX13, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,PEX13Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, PEX13Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,PEX13gain)

PIN3 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "PIN3 in ABP1")
PIN3Lost = anti_join(ABP1inABP1,PIN3, by = "DHFR3")
PIN3Lost %<>% ungroup()
PIN3Lost%<>% mutate(DHFR1.2 = "PIN3 in ABP1")
PIN3Lost%<>% group_by(DHFR1.2)
PIN3Ratio = anti_join(ABP1inABP1,PIN3Lost, by = "DHFR3")
PIN3Ratio = merge(PIN3,PIN3Ratio, by = "DHFR3")
PIN3Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
PIN3Ratio[,c(6:9)]=NULL
colnames(PIN3Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
PIN3gain = anti_join(PIN3, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,PIN3Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, PIN3Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,PIN3gain)

RVS167 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "RVS167 in ABP1")
RVS167Lost = anti_join(ABP1inABP1,RVS167, by = "DHFR3")
RVS167Lost %<>% ungroup()
RVS167Lost%<>% mutate(DHFR1.2 = "RVS167 in ABP1")
RVS167Lost%<>% group_by(DHFR1.2)
RVS167Ratio = anti_join(ABP1inABP1,RVS167Lost, by = "DHFR3")
RVS167Ratio = merge(RVS167,RVS167Ratio, by = "DHFR3")
RVS167Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
RVS167Ratio[,c(6:9)]=NULL
colnames(RVS167Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
RVS167gain = anti_join(RVS167, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,RVS167Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, RVS167Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,RVS167gain)

SDC25 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "SDC25_17W in ABP1")
SDC25Lost = anti_join(ABP1inABP1,SDC25, by = "DHFR3")
SDC25Lost %<>% ungroup()
SDC25Lost%<>% mutate(DHFR1.2 = "SDC25_17W in ABP1")
SDC25Lost%<>% group_by(DHFR1.2)
SDC25Ratio = anti_join(ABP1inABP1,SDC25Lost, by = "DHFR3")
SDC25Ratio = merge(SDC25,SDC25Ratio, by = "DHFR3")
SDC25Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
SDC25Ratio[,c(6:9)]=NULL
colnames(SDC25Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
SDC25gain = anti_join(SDC25, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,SDC25Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, SDC25Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,SDC25gain)

SHO1 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "SHO1 in ABP1")
SHO1Lost = anti_join(ABP1inABP1,SHO1, by = "DHFR3")
SHO1Lost %<>% ungroup()
SHO1Lost%<>% mutate(DHFR1.2 = "SHO1 in ABP1")
SHO1Lost%<>% group_by(DHFR1.2)
SHO1Ratio = anti_join(ABP1inABP1,SHO1Lost, by = "DHFR3")
SHO1Ratio = merge(SHO1,SHO1Ratio, by = "DHFR3")
SHO1Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
SHO1Ratio[,c(6:9)]=NULL
colnames(SHO1Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
SHO1gain = anti_join(SHO1, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,SHO1Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, SHO1Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,SHO1gain)

SLA11 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "SLA1-1 in ABP1")
SLA11Lost = anti_join(ABP1inABP1,SLA11, by = "DHFR3")
SLA11Lost %<>% ungroup()
SLA11Lost%<>% mutate(DHFR1.2 = "SLA1-1 in ABP1")
SLA11Lost%<>% group_by(DHFR1.2)
SLA11Ratio = anti_join(ABP1inABP1,SLA11Lost, by = "DHFR3")
SLA11Ratio = merge(SLA11,SLA11Ratio, by = "DHFR3")
SLA11Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
SLA11Ratio[,c(6:9)]=NULL
colnames(SLA11Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
SLA11gain = anti_join(SLA11, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,SLA11Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, SLA11Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,SLA11gain)

SLA12 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "SLA1-2 in ABP1")
SLA12Lost = anti_join(ABP1inABP1,SLA12, by = "DHFR3")
SLA12Lost %<>% ungroup()
SLA12Lost%<>% mutate(DHFR1.2 = "SLA1-2 in ABP1")
SLA12Lost%<>% group_by(DHFR1.2)
SLA12Ratio = anti_join(ABP1inABP1,SLA12Lost, by = "DHFR3")
SLA12Ratio = merge(SLA12,SLA12Ratio, by = "DHFR3")
SLA12Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
SLA12Ratio[,c(6:9)]=NULL
colnames(SLA12Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
SLA12gain = anti_join(SLA12, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,SLA12Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, SLA12Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,SLA12gain)

SLA13 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "SLA1-3 in ABP1")
SLA13Lost = anti_join(ABP1inABP1,SLA13, by = "DHFR3")
SLA13Lost %<>% ungroup()
SLA13Lost%<>% mutate(DHFR1.2 = "SLA1-3 in ABP1")
SLA13Lost%<>% group_by(DHFR1.2)
SLA13Ratio = anti_join(ABP1inABP1,SLA13Lost, by = "DHFR3")
SLA13Ratio = merge(SLA13,SLA13Ratio, by = "DHFR3")
SLA13Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
SLA13Ratio[,c(6:9)]=NULL
colnames(SLA13Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
SLA13gain = anti_join(SLA13, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,SLA13Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, SLA13Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,SLA13gain)

YSC84 = subset(ABP1Swap,ABP1Swap$DHFR1.2 == "YSC84 in ABP1")
YSC84Lost = anti_join(ABP1inABP1,YSC84, by = "DHFR3")
YSC84Lost %<>% ungroup()
YSC84Lost%<>% mutate(DHFR1.2 = "YSC84 in ABP1")
YSC84Lost%<>% group_by(DHFR1.2)
YSC84Ratio = anti_join(ABP1inABP1,YSC84Lost, by = "DHFR3")
YSC84Ratio = merge(YSC84,YSC84Ratio, by = "DHFR3")
YSC84Ratio %<>% mutate(Ratio = medrep.x/medrep.y)
YSC84Ratio[,c(6:9)]=NULL
colnames(YSC84Ratio) = c("DHFR3","DHFR1.2","n()","medsize","medrep","Ratio")
YSC84gain = anti_join(YSC84, ABP1inABP1, by ="DHFR3")
ABP1SWAPLost = rbind(ABP1SWAPLost,YSC84Lost)
ABP1SWAPRatio = rbind(ABP1SWAPRatio, YSC84Ratio)
ABP1SWAPGain = rbind(ABP1SWAPGain,YSC84gain)

# The threshold for PPIs that are stronger or weaker following the domain swapping correspond to the top and bottom 7.3%.

quantile(ABP1SWAPRatio$Ratio, .927)
quantile(ABP1SWAPRatio$Ratio, .073)

# Top 7.3% is 1.319 and bot 7.3% is 0.639.
# Distribution of Ratios with the thresholds.

DistributionRatioABP1SWAP = ggplot(data = ABP1SWAPRatio, aes(x = Ratio, colour = DHFR1.2)) + scale_x_continuous(breaks=seq(-10,7,0.5)) + geom_density() + geom_vline(xintercept=c(0,0.639,1.319))
pdf("DistributionRatioABP1SWAP_0120Analysis.pdf")
DistributionRatioABP1SWAP + theme(legend.position = "none")
dev.off()

ABP1SWAPRatioPositive = subset(ABP1SWAPRatio,ABP1SWAPRatio$Ratio >= 1.319)
ABP1SWAPRatioPositive%<>% group_by(DHFR1.2)
ABP1SWAPRatioNegative = subset(ABP1SWAPRatio, ABP1SWAPRatio$Ratio <= 0.639)
ABP1SWAPRatioNegative%<>% group_by(DHFR1.2)

# Combine everything in one dataframe. Lost PPIs are assigned a Ratio 0 and gained ones are assigned a Ratio of 10.

ABP1SWAPLost%<>% mutate(Ratio = 0)
ABP1SWAPLost%<>% group_by(DHFR1.2)
ABP1SWAPGain%<>% mutate(Ratio = 10)
ABP1SWAPGain%<>% group_by(DHFR1.2)

ABP1SWAPIntNegAffected = rbind(ABP1SWAPLost,ABP1SWAPRatioNegative)
ABP1SWAPIntPosAffected = rbind(ABP1SWAPGain,ABP1SWAPRatioPositive)

ABP1SWAPIntAffected = rbind(ABP1SWAPIntNegAffected,ABP1SWAPIntPosAffected)
write.csv(ABP1SWAPIntAffected,"PCAOct18ABP1SWAPALLIntAffec_0120Analysis.csv")



