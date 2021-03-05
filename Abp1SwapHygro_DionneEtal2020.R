library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(growthcurver)
library(ggrepel)
library(forcats)

# Script for the analysis of the growth of the different Abp1 SH3-swapped strains in liquid media containing Hygromycin B.

# The raw data available with the article correspond to this step. Some colunm names and some strain names have been edited manually for the paper.

dat <- read_tsv("SH3paperData - FigS3G-H.Abp1SwapHygroGrowthcurves.tsv")
table(dat$strain)
dat <- dat %>% rename(well = well.positions) 

# Data comes from two plates, add plate id.

dat$plate <- c(rep(1,96),rep(2,nrow(dat)-96))

# Change format of the dataframe.

datl <- dat %>% mutate_if(is.numeric,as.character, is.factor, as.character) %>% 
  pivot_longer(cols=4:195, names_to="time", values_to="od")

datl$od  <- as.numeric(as.character(datl$od))
datl$time  <- as.numeric(as.character(datl$time))*60 #change to minutes

# Remove datapoint above 1.1.

datl %<>% mutate(od = ifelse(od>1.1,NA,od)) %>% filter(!is.na(od)) 

# Calculate the growth rate.

growth_rates <- datl %>% group_by(condition, strain, well,plate) %>%
  mutate(note_fit = SummarizeGrowth(time, od)$vals$note,
         rval = SummarizeGrowth(time, od)$vals$r,
         kval = SummarizeGrowth(time, od)$vals$k,
         tmid = SummarizeGrowth(time, od)$vals$t_mid,
         aucexp = SummarizeGrowth(time, od)$vals$auc_l,
         max_size = nth(od, -3)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup()


# Add information on similarity of sequences (Identity) of yeast SH3 compared to Abp1 SH3. This data is available with the article.

dis <- read.table("SH3_seq_indentity_n31.txt", header=T, sep=" ",row.names = NULL)
sim_to_abp1 <- data.frame(dis$row.names, dis$Abp1)
names(sim_to_abp1) <- c("dom", "sim")
sim_to_abp1$dom <- toupper(sim_to_abp1$dom)

# Merge growth data with sequence identity.

com_data <- left_join(growth_rates,sim_to_abp1, by=c("strain"="dom"))

# Annotate some strains.

humans = c("CTTN", "HCLS1", "DBNL")
controls = c("ABP1", "ABP1 delta SH3", "ABP1 WT", "BY4741","delta ABP1")

com_data %<>% mutate(sim_mod = ifelse(strain=="ABP1 delta SH3", 0,
                                      ifelse(strain=="BY4741", 1,
                                             ifelse(strain=="ABP1 WT", 1,
                                                    ifelse(strain=="delta ABP1", 0, sim)))),
                     origin_SH3 = ifelse(strain %in% humans,"human ABP1 ortholog",
                                         ifelse(strain %in% controls, "control","yeast SH3")))


# Keep only the hygromycin B condition.

hygro <- com_data %>% filter(condition=="Hygro", !is.na(sim_mod))
write.csv(hygro,"Abp1SWAPHygroCRLAnalysis_0320.csv")


hygro2 =read.csv ("Abp1SWAPHygroCRLAnalysis_0320.csv")
hygro2$condition <- (as.character(hygro2$condition))
hygro2$strain <- (as.character(hygro2$strain))
hygro2$origin_SH3 <- (as.character(hygro2$origin_SH3))


# Boxplot corresponding to Figure S2E. The figure was also edited in illustrator.

pdf("BoxPlotAbp1SwapHygro_CRLCode_0320.pdf",7, 10)
ggplot(hygro2, aes(x=reorder(strain, 1-sim_mod, fun = median),y=aucexp, fill=origin_SH3)) +
  geom_boxplot(outlier.size=-1) +
  theme_bw()+
  scale_fill_manual(values=c("grey1","grey39","grey78"))+
  labs(y = "Growth in hygromycin (AU)", x="Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20))+
  coord_flip()+
  theme(legend.position = c(0.2, 0.8))

dev.off()
