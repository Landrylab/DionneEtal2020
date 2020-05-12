library(dplyr)
library(rowr)
library(magrittr)
library(gplots)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(stringr)
library(UpSetR)

# Script for the analysis of the liquid DHFR-PCA experiment for the validation of the quantitative changes detected upon SH3-deletions in the solid DHFR-PCA screen.
# Get 96-well plate definitions and optic density data. _2 is for PPIs that were originally detected as gained and _1 is for lost.

Plate1qA= read.csv("plate1_qA.csv", header = T, sep= ",")
MAPplate1qA = read.csv("plate1_qA_position.csv", header = T)

Plate1qB = read.csv("plate1_qB.csv", header =T, sep= ",")
MAPplate1qB = read.csv("plate1_qB_position.csv", header = T)

Plate2qA = read.csv("plate2_qA.csv", header =T, sep=",")
MAPplate2qA = read.csv("plate2_qA_position.csv", header = T)

Plate2qB = read.csv("plate2_qB.csv", header = T, sep=",")
MAPplate2qB = read.csv("plate2_qB_position.csv", header = T)

Plate3qA = read.csv("plate3_qA.csv", header = T, sep=",")
MAPplate3qA = read.csv("plate3_qA_position.csv", header = T)

Plate3qB = read.csv("plate3_qB.csv", header = T, sep=",")
MAPplate3qB = read.csv("plate3_qB_position.csv", header = T)

Plate4qA = read.csv("plate4_qA.csv", header = T, sep=",")
MAPplate4qA = read.csv("plate4_qA_position.csv", header = T)

Plate4qB = read.csv("plate4_qB.csv", header = T, sep=",")
MAPplate4qB = read.csv("plate4_qB_position.csv", header = T)

# Merge 96-well plate definitions and optic density data.

Plate1qAMerge = merge(Plate1qA,MAPplate1qA, by = "Well.positions")
Plate1qAMerge = subset(Plate1qAMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate1qAMerge[,1]=NULL

Plate1qBMerge = merge(Plate1qB,MAPplate1qB, by = "Well.positions")
Plate1qBMerge = subset(Plate1qBMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate1qBMerge[,1]=NULL

Plate2qAMerge = merge(Plate2qA,MAPplate2qA, by = "Well.positions")
Plate2qAMerge = subset(Plate2qAMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate2qAMerge[,1]=NULL

Plate2qBMerge = merge(Plate2qB,MAPplate2qB, by = "Well.positions")
Plate2qBMerge = subset(Plate2qBMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate2qBMerge[,1]=NULL

Plate3qAMerge = merge(Plate3qA,MAPplate3qA, by = "Well.positions")
Plate3qAMerge = subset(Plate3qAMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate3qAMerge[,1]=NULL

Plate3qBMerge = merge(Plate3qB,MAPplate3qB, by = "Well.positions")
Plate3qBMerge = subset(Plate3qBMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate3qBMerge[,1]=NULL

Plate4qAMerge = merge(Plate4qA,MAPplate4qA, by = "Well.positions")
Plate4qAMerge = subset(Plate4qAMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate4qAMerge[,1]=NULL

Plate4qBMerge = merge(Plate4qB,MAPplate4qB, by = "Well.positions")
Plate4qBMerge = subset(Plate4qBMerge, select=c(Well.positions,Bait,Prey,X0s:X258333s))
Plate4qBMerge[,1]=NULL

AllData = rbind(Plate1qAMerge, Plate1qBMerge, Plate2qAMerge, Plate2qBMerge, Plate3qAMerge, Plate3qBMerge, Plate4qAMerge, Plate4qBMerge)
colnames(AllData) = c("Bait","Prey",seq(0,287,1))
write.csv(AllData,"LiquidPCAAllDATA_0320.csv")
# This corresponds to the raw dataset available with the article. Some colunm names and some bait names have been edited manually for the paper.

# Adapt the OD to 1 mL of culutre.

AllDataNorm = AllData * 4
AllDataNorm = cbind(AllData[,1:3], AllDataNorm[,4:290])
Empty = AllDataNorm[grep("Empty", AllDataNorm$Bait),]
AllDataFinal = anti_join(AllDataNorm,Empty)

# comparison of the last time point of WT Vs stuffed interactions

LastTP = subset(AllDataFinal, select=c(Bait,Prey,290))
Stuffed = LastTP[grep("Stuffed", LastTP$Bait),]
Stuffed%<>% mutate(SH3 = "Stuffed")
WT = anti_join(LastTP,Stuffed)
WT%<>% mutate(SH3 = "WT")

setDT(Stuffed)[, paste0("Bait", 1:2) := tstrsplit(Bait, " ")]
Stuffed = subset(Stuffed, select=c(Bait1,Prey,3,SH3))
colnames(Stuffed) = c("Bait","Prey","287","SH3")

WT$Prey <- (as.character(WT$Prey))
Stuffed$Prey <- (as.character(Stuffed$Prey))

WTVsStuffed = merge(WT,Stuffed, by = c("Bait", "Prey"))
WTVsStuffed%<>% mutate(Ratio = `287.y`/`287.x`)

WTVsStuffed = subset(WTVsStuffed, select = c(Bait,Prey,`287.x`,`287.y`,Ratio))
colnames(WTVsStuffed) = c("Bait","Prey","WTscore","Stuffedscore","Ratio")

Gained = WTVsStuffed[grep("_2", WTVsStuffed$Prey),]
Gained%<>% mutate(InteractionType = "gained")

Lost = WTVsStuffed[grep("_1", WTVsStuffed$Prey),]
Lost%<>% mutate(InteractionType = "Lost")

WTVsStuffedAnnotated = rbind(Gained, Lost)
write.csv(WTVsStuffedAnnotated,"PCAliquidValidation_WTvsStuffedAnn_0220Analysis.csv")
# This corresponds to the analyzed dataset available with the article. Some colunm names and some bait names have been edited manually for the paper.

## Figure S2C.
#boxplot with last timepoint OD.

PCAliquidValid = read.csv("PCAliquidValidation_WTvsStuffedAnn_0220Analysis.csv")

my_paletteS2 = colorRampPalette(c("#B2182B","#2166AC"))(n = 2)

pdf("BoxPLotPCAliquid_0220Analysis_Fig0320.pdf",10, 10)
p = boxplot(Ratio~InteractionType, data = PCAliquidValid, main = "PCA validation", ylab = "OD Ratios Stuffed/WT", par(cex.axis=1.5), las = 1,outline = FALSE, col = my_paletteS2,ylim = c(0, 2.5),par(cex.lab=2.0), boxwex = 0.6, par(mar=c(5.0,6.0,4.0,2.0)))
text(x=c(1:2), y=p$stats[nrow(p$stats),] + 0.1, paste("n = ",table(PCAliquidValid$InteractionType),sep=""), cex = 2.0, col = "black") 
abline(h=1.0, col="black",lty=5)
dev.off()
UniqueGained = subset(PCAliquidValid,PCAliquidValid$InteractionType == "gained")
UniqueLost = subset(PCAliquidValid,PCAliquidValid$InteractionType == "Lost")


## Figure S2D
#Growth curves of Pbs2/YJL128C interaction. 

Pbs2Data = AllDataFinal[grep("YJL128C",AllDataFinal$Prey),]

colours = c("SHO1" = "#F26C4F", "SHO1 Stuffed" = "#FFF568", "NBP2" = "#ACD373", "NBP2 Stuffed" = "#00BFF3")
reshape.subdf = melt(Pbs2Data)
names(reshape.subdf) = c("bait", "prey", "time", "OD")
pdf("GrowthCurvesPBS2PPI_0420.pdf")
reshape.subdf %>% ggplot(aes(time, OD, colour = bait)) + theme_bw() + ylab("O.D. 600 nm") + scale_x_discrete("Time (mins)", breaks=seq(0,287,50)) + geom_point(size = 1) + scale_colour_manual(values = colours)
dev.off()
