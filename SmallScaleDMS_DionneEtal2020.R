library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(growthcurver)
library(ggrepel)
library(forcats)

# Script for the analysis of the small scale DMS experiment (selected Abp1 SH3 mutants, Abp1-Hua2 DHFR-PCA) and its comparison with the large scale DMS results.

dat <- read_tsv("SH3paperData - Fig4BCD.PCAgrowhtSelectedDMSmutants.tsv")
dat %<>% rename(strain = Strain) 

# Add sample id.

dat$id <- 1:nrow(dat)

# Change the format of the dataframe.

datl <- dat %>% mutate_if(is.numeric,as.character, is.factor, as.character) %>% 
  pivot_longer(cols=2:289, names_to="time", values_to="od")

datl$od  <- as.numeric(as.character(datl$od))
datl$time  <- as.numeric(as.character(datl$time))*15 #change to minutes

# This dataset correspond to the raw data that is available with the article. Some column names and strain names have been changed for the paper.
write.csv(datl,"DMSValidationRaw_0420.csv")

# Get the DMS data. The data used correspond to the ratio..log2.scale..mean from the DMS data normalized to the Reference condition.
# This data is available with the article. Some column names and strain names have been changed for the paper.

dms <- read.table("SH3paperData - Fig4.DMSData.tsv", sep="\t", header=T)

# Keep only a.a and HUA2_MTX2.

dms %<>% filter(condition=="HUA2_MTX2", mutation.format=="amino acid")

# Separate the position.reference column.

dms %<>% separate(position.reference,c("pos", "wtaa"), sep=" ") %>%
  mutate(pos=as.integer(pos)+1)

# Calculate growth in small scale experiment (AUC).

growth <- datl %>% group_by(strain, id) %>%
  mutate(aucexp = sum(od)) %>%
  slice(1) %>%
  ungroup()

# Calculate the average between replicates.

growth_ave <- growth %>% group_by(strain) %>%
  summarize(ave_au = mean(aucexp))

# Merge with DMS data. Create ID equivalent to strain ID to merge the two types of data.

dms %<>% mutate(id_strain = paste(wtaa,pos,mutated, sep=""))

com<- left_join(growth_ave,dms, by=c("strain"="id_strain"))

# Figure S4D. The figure was also edited in illustrator.

ref_wt <- com$ave_au[com$strain=="WT"]
ref_del <- com$ave_au[com$strain=="stuffed"]

test_cor <- cor.test(com$ave_au,com$ratio..log2.scale..mean, method="spearman")

correlation_on_plot = paste("rho=", round(test_cor$estimate[1],2), " ", "p-value = 1e-13")                           

pdf("DMSvalidSelecMutPCA_CRLCode_0320.pdf",12, 10)
ggplot(com, aes(x=ratio..log2.scale..mean, y=ave_au))+
  geom_point()+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),) +
  labs(y = "Growth in MTX", x="DMS score") +
  theme_bw()+
  geom_hline(yintercept=ref_del, linetype="dashed", 
             color = "red", size=1)+
  geom_hline(yintercept=ref_wt, linetype="dashed", 
             color = "blue", size=1)+
  geom_label_repel(aes(label = strain))+
  annotate("text", x=-1, y=ref_wt+3, label= "WT", col="blue", size=7)+
  annotate("text", x=3, y=ref_del-3, label= "Deletion", col="red", size=7)+
  annotate("text", x= 0, y= 110, label = correlation_on_plot, size=5)

dev.off()

# Dataset of analyzed data available with the paper correspond to this dataframe. Some column names and bait names have been manually changed for the paper.
write.csv(com,"DMSValidationAnalyzed_0420.csv")