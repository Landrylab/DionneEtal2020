source("http://bioconductor.org/biocLite.R")
biocLite("EBImage")
install.packages('gitter')
require(gitter)
library(gitter)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(forcats)
library(directlabels)
library(heatmaply)
library(RColorBrewer)
library(cowplot)
library(gplots)
library(plotly)

# Script for the analysis of the growth screen on solid media in stress conditions.

# Gitter analysis with SC plate 1 last timepoint as REF.

Ref="C:/Users/ugols/Desktop/GrowthScreen/37degre/d000117_370_001_19-04-14_15-14-52.JPG"
f="C:/Users/ugols/Desktop/GrowthScreen/37degre"
path_dat="C:/Users/ugols/Desktop/GrowthScreen/37degre/37degre_gitter"
path_grid="C:/Users/ugols/Desktop/GrowthScreen/37degre/37degre_gitter"
gitter.batch(image.files=f,ref.image.file=Ref,plate.format=1536, verbose='p',grid.save=path_grid, dat.save=path_dat)

# Merge Plate definitions, strain positions on array with the Gitter analysis.

PlatesPosition = read.csv("Plate_Positions_170419.csv", header = T, sep = ";")
platesDefinition37Deg = read.csv("Plate_Definition_37degres_240419.csv", header = T, sep = ";")
GS37 = data.frame(matrix(vector(), 0, 10, dimnames = list(c(), c("row", "col", "size", "circularity", "flags", "orf", "Plate", "Condition", "Time", "Rep"))))
for (i in 1:nrow(platesDefinition37Deg)) {
  data = read.table(paste("37degre_gitter/",as.character(platesDefinition37Deg[i,"dat_Gitter"]), sep=""), sep="\t", header=F, stringsAsFactor=F, strip.white = T,comment.char = "#", na.strings="")
  names(data) <- c("row", "col", "size", "circularity", "flags")
  Plate = platesDefinition37Deg[i, "Plate"]
  Condition = platesDefinition37Deg[i, "Condition"]
  Time = platesDefinition37Deg[i, "Time"]
  Rep = platesDefinition37Deg[i, "Rep"]
  subdf <- cbind(data, rep(Plate, nrow(data)), rep(Condition, nrow(data)), rep(Time, nrow(data)), rep(Rep, nrow(data)))
  subdf <- merge(subdf, PlatesPosition[,c("orf", "row", "col")], by.x=c("col", "row"), by.y=c("col", "row"))
  names(subdf) <- c("row", "col", "size", "circularity", "flags","orf", "Plate", "Condition", "Time", "Rep")
  GS37 = rbind(GS37, subdf)
}
write.csv(GS37,"AlldataGS37_280619.csv")

colnames(GS37) = c("col", "row", "size", "circularity","flags","Plate","Condition","Time","Rep","orf")
GS37$orf = as.character(GS37$orf)
write.csv(GS37,"AlldataGS37_ColOK_280619.csv")
# The dataset of raw data available with the article correspond to this step. Some colunm names and some strain names have been edited manually for the paper.

dat <- read_delim("AlldataGS37_ColOK_280619.csv", delim=",")

# Log 2 transformation of the data.

dat %<>% mutate(log2 = log2(size+0.5))

# Remove abnormal colonies by using the flags called by Gitter based on the circularity of the colonies. Also, remove the borders.

dat %<>% filter(row!=1, row!=2, row!=32,row!=31,  
                col!=1, col!=2, col!=47,col!=48) %>%
  mutate(flags = ifelse(is.na(flags),"OK", flags))  %>% 
  filter(flags !="C",flags !="S,C")

# Calculate median per plate per time point to determine the timepoint to use for the analysis and correct for plate effects.

fin_size <- dat %>% group_by(Condition, Rep, Time) %>%
  summarize(med_size = median(log2))

# Distribution of median size per plate through time to choose time to consider.

fin_size %>%  
  ggplot(., aes(x=Time, y=med_size, col=as.factor(interaction(Condition,Rep)))) +
  geom_line()+
  theme(legend.position = "none")+
  geom_dl(aes(label = interaction(Condition,Rep)), method = list(dl.combine("first.points", "last.points"), cex = 0.8))

# Remove conditions with no growth or with heterogeneity of drug distribution in the plate media.

bad_conditions <- c("Nitroquinoline 200uM","Rapamycin 150nM","Nystatin 1uM","Caffein 30mM","CuSO4 10mM",
                    "NaCl 2.5M")
dat %<>% filter(!(Condition %in% bad_conditions))

# Use time = 74 because most conditions stop growthing after that.

times_to_keep = c(0,74)
datf <- dat %>% filter( Time %in% times_to_keep) %>% select(row, col, Condition,Plate, orf, Time, Rep, log2) %>%
  pivot_wider(names_from=c(Time), values_from = log2) %>%
  rename(tini = "0", tfin="74")

# Correcting for the difference between replicated plates will correct for the slightly different initial conditions.
# Put replicates side by side.

datw <- datf %>% select(row, col, Condition, orf, tfin, Rep) %>%
  pivot_wider(names_from=c(Rep), values_from = tfin) %>%
  rename(rep1= "1", rep2="2")

# Comparing plate replicates.

datw %>% ggplot(., aes(x=rep1, y=rep2, color=Condition))+
  geom_point()

cor.test(datw$rep1, datw$rep2)

# Calculate median size per condition to correct the two replicates.
# Correct plate replicate 2 so its median value per condition is the same as plate replicate 1.

datw %<>% group_by(Condition) %>%
  mutate(med_per_cond_rep_tfin1 = median(rep1, na.rm=T),
         med_per_cond_rep_tfin2 = median(rep2, na.rm=T)) %>%
  ungroup() %>%
  mutate(corrected_rep2 = rep2 + (med_per_cond_rep_tfin1 - med_per_cond_rep_tfin2))

# Take median among replicates in a given plate replicate and remove points where replicates disagree more than 2-fold (data in log2) between plate replicates.

datw_ave <- datw %>% group_by(Condition, orf) %>%
  summarize(med_orf_per_rep1 = median(rep1, na.rm=T),
            med_orf_per_rep2 = median(corrected_rep2, na.rm=T)) %>%
  filter(abs(med_orf_per_rep1-med_orf_per_rep2) < 2)


# Calculate average of the two plate replicates (which correspond to medians of within plate values).

datw_ave %<>% rowwise() %>%
  mutate(mean_score = mean(c(med_orf_per_rep1,med_orf_per_rep2), na.rm=T))



## Analysis of the growth of Abp1 SH3-swapped strains.

list_orfs <- datw_ave %>% dplyr::filter(grepl('ABP1|Abp1|abp1', orf)) %>% 
  select(orf)
list_orfs <- unique(list_orfs$orf)

# Removing Abp1 SH3 in Sla1 strains.

list_orfs <- list_orfs [c(-34,-35,-36)]

# Add information on similarity of sequences (Identity) of yeast SH3 compared to Abp1 SH3. This data is available with the article.

dis <- read.table("SH3_seq_indentity_n31.txt", header=T, sep=" ",row.names = NULL)
sim_to_abp1 <- data.frame(dis$row.names, dis$Abp1)
names(sim_to_abp1) <- c("dom", "sim")
sim_to_abp1$dom <- toupper(sim_to_abp1$dom)

# Create a matrix to add the sequence similarity scores.

mat <- datw_ave %>% filter(orf %in% list_orfs) %>%
  select(orf, Condition, mean_score) %>%
  pivot_wider(names_from=c(Condition), values_from = mean_score)

# Modification of strain names.

mat %<>% mutate(orf = ifelse(orf=="ABP1", "ABP1 WT",
                             ifelse(orf=="ABP1 Stuffed", "ABP1 delta SH3",
                                    ifelse(orf=="ABP1 KO", "ABP1 delta",orf))))

mat %<>% mutate (orf = str_replace(orf, " in ABP1",""),
                 orf = str_replace(orf, "SDC25_17W","SDC25"))

# Merge growth and sequence identity of SH3.

mat_dist <- left_join(mat, sim_to_abp1, by=c("orf"="dom"))

# Annotation of strains not SH3 swapped.

mat_dist %<>% mutate(sim = ifelse(orf=="ABP1 WT", 1,
                                  ifelse(orf == "ABP1 delta", 0,
                                         ifelse (orf == "ABP1 delta SH3",0,sim))))

# Create bins of levels of conservation for the Heatmap.

mat_dist %<>% mutate(Conservation_class = ifelse(sim==0, "No domain",
                                                 ifelse(sim >= 0.1 & sim < 0.2,"[0.1,0.2[",
                                                        ifelse(sim >= 0.2 & sim < 0.3,"[0.2,0.3[",
                                                               ifelse(sim >= 0.3 & sim < 0.4,"[0.3,0.4[",
                                                                      ifelse(sim >= 0.4 & sim < 0.5,"[0.4,0.5[",
                                                                             ifelse(sim >= 0.5 & sim < 0.6,"[0.5,0.6[",
                                                                                    ifelse(sim == 1,"ABP1 SH3","N"))))))))       

mat_dist <- data.frame(mat_dist, row.names = mat_dist$orf)

# This correspond to the analyzed data for the Abp1 SH3-swapped growth that is available with the article. Some colunm names and some strain names have been edited manually for the paper. 

#Heatmap

seq_sim <- mat_dist$Conservation_class

mat_dist <- mat_dist[,2:(ncol(mat_dist)-2)]

pal <- brewer.pal(n = 8, name = "Greys")

pal_col <- c("No domain"= pal[1], "[0.1,0.2[" = pal[2],"[0.2,0.3[" = pal[3],
             "[0.3,0.4[" = pal[4],
             "[0.4,0.5["= pal[5],
             "[0.5,0.6["= pal[6],
             "ABP1 SH3" = pal[7])

# Heatmap of data, scaling the rows. If we would assume all variables come from some normal distribution,  
# then scaling (i.e.: subtract the mean and divide by the standard deviation) would bring them all close to the standard normal distribution. 
# In such a case, each value would reflect the distance from the mean in units of standard deviation.
# Figure S3E. The figure was also edited in illustrator.

dir.create("folder")
heatmaply_cor(
  mat_dist,
  xlab = "Condition",
  ylab = "Genotype",
  scale = "row",
  k_col = 1,
  k_row = 1,
  plot_method = "plotly",
  RowSideColors = seq_sim,
  row_dend_left = TRUE,
  row_side_palette= pal_col,
  side_color_colorbar_len = 0.2, file = "folder/Abp1SWAPHygroAndSim_CRLcode_0320.html")


#### Comparison of the growth of the SH3-deleted or KO strains with the WT strains.
# Use a table of definitions of the strains for the annotation. This table is available with the article.

geno <- read_delim("NameConditionsSolidKO_D_WT_VUD.txt", delim="\t")
datw_ave <- left_join(datw_ave,geno, by=c("orf"="orf"))

# select stuffed, WT and KO strains.
# Remove YSC84 because we don't have the SH3-deleted (stuffed) strain.
# Remove CDC25 because there is no KO strain as it is the only SH3-containing gene that is essential.
# For multi-SH3 proteins, only keep the WT, KO and all SH3-deleted strains.

sel_construct = c("KO", "D", "WT")
remove = c("BEM1-1 Stuffed","BEM1-2 Stuffed","BZZ1-1 Stuffed", "BZZ1-2 Stuffed",
           "SLA1 -|-|WT","SLA1 -|WT|-","SLA1 -|WT|WT","SLA1 WT|-|WT","SLA1 WT|WT|-",
           "SLA1 WT|-|-", "YSC84 KO","CDC25")
           
comp_ko_d <- datw_ave %>% filter(type_construct %in% sel_construct) %>%
  filter(!(orf %in% remove))

# Some names in gene_tagged have space and cause problems.

comp_ko_d %<>% mutate(gene_tagged = str_replace_all(gene_tagged, "[^[:alnum:]]", ""))

# Calculate difference of growth between SH3-deleted (stuffed) or KO strains with the WT strains. Do this per condition.

diff_with_WT <- comp_ko_d %>% select(gene_tagged, Condition, mean_score,type_construct) %>%
  pivot_wider(names_from=c(type_construct), values_from = mean_score) %>%
  rowwise() %>%
  mutate(KO_minus_WT = KO-WT, D_minus_WT = D-WT) %>%
  select(gene_tagged,Condition,WT_minus_KO,WT_minus_D) %>%
  pivot_longer(cols=3:4, names_to="type_diff", values_to="diff")

# Boxplots.

f2d <- diff_with_WT %>% filter(!is.na(diff),gene_tagged !="HOF1")
pdf("BoxPlotStuffedKOVsWT_CRLCode_0320.pdf",7, 10)
f2d%>%
  ggplot(., aes(x=reorder(gene_tagged, diff,fun=median), y=diff, fill=type_diff)) +
  geom_boxplot(outlier.size=-1)+ 
  theme_bw()+
  ylim(-5,1)+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'))+
  scale_fill_manual(values=c("grey1","grey78"),name = "Comparison versus WT", labels = c("Deletion of SH3","Gene deletion"))+
  labs(y="Growth relative to WT", x="Gene")+
  coord_flip()+
  theme(legend.position = c(0.2, 0.8))
dev.off()

# Analyze growth of strains with multi-SH3 containing proteins.

conserved = c("BEM1-1 Stuffed",
              "BEM1-2 Stuffed",
              "BEM1 KO",
              "BEM1 doubled Stuffed",
              "BEM1",
              "BZZ1 KO",
              "BZZ1-1 Stuffed", 
              "BZZ1-2 Stuffed",
              "BZZ1",
              "BZZ1 doubled Stuffed",
              "SLA1 triple Stuffed",
              "SLA1 WT|WT|-",
              "SLA1 WT|-|WT",
              "SLA1 -|WT|WT",
              "SLA1 WT|WT|WT",
              "SLA1 KO")

# Change Sla1 nomenclature.

mSH3_ko_d <- datw_ave %>% filter(orf %in% conserved) %>%
  mutate(orf = str_replace(orf, "SLA1 WT\\|WT\\|WT","SLA1"))


# Growth of WT strains per condition.

mSH3_ko_d_WTscores <-  mSH3_ko_d %>% filter(type_construct=="WT") %>% 
  select(Condition, gene_tagged, mean_score) %>%
  rename(WT = mean_score)

# Growth of SH3-deleted and KO strains per condition.

mSH3_ko_d_scores   <-  mSH3_ko_d %>% filter(type_construct !="WT")                      

# Merge the two.

mSH3_ko_d_scores <- left_join(mSH3_ko_d_scores,mSH3_ko_d_WTscores, 
                              by=c("Condition"="Condition",
                                   "gene_tagged"="gene_tagged")) %>%
  mutate(rel_to_wt = mean_score-WT)

# Add a rank  to be able to order them in the plot.

mSH3_ko_d_scores %<>% mutate( rank = 
                                ifelse(orf == "BZZ1 KO",1,
                                       ifelse(orf == "BZZ1-1 Stuffed",2,
                                              ifelse(orf == "BZZ1-2 Stuffed",3,
                                                     ifelse(orf == "BZZ1 doubled Stuffed",4,
                                                            ifelse(orf == "BEM1 KO",5,
                                                                   ifelse(orf == "BEM1-1 Stuffed",6,       
                                                                          ifelse(orf == "BEM1-2 Stuffed",7,
                                                                                 ifelse(orf == "BEM1 doubled Stuffed",8,
                                                                                        ifelse(orf == "SLA1 KO",9,
                                                                                               ifelse(orf == "SLA1 -|WT|WT",10,
                                                                                                      ifelse(orf == "SLA1 WT|-|WT",11,
                                                                                                             ifelse(orf == "SLA1 WT|WT|-",12,
                                                                                                                    ifelse(orf == "SLA1 triple Stuffed",13,
                                                                                                                           0)))))))))))))) 

# Boxplots.

pdf("BoxPlotStuffedKOVsWTMSH3_CRLCode_0320.pdf",7, 10)
ggplot(mSH3_ko_d_scores, aes(x=reorder(orf, desc(rank),fun=median), y=rel_to_wt, fill=type_construct))+
  geom_boxplot(outlier.size=-1)+
  theme_bw()+
  ylim(-5,1)+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'))+
  scale_fill_manual(values=c("grey1","grey78"),name = "Comparison versus WT", 
                    labels = c("Deletion of SH3","Gene deletion"))+
  labs(y="Growth relative to WT", x="Gene")+
  scale_x_discrete(labels = c("D|D|D","WT|WT|D","WT|D|WT","D|WT|WT","SLA1","D|D","WT|D","D|WT","BEM1","D|D","WT|D","D|WT","BZZ1") )+
  coord_flip()+
  theme(legend.position = "none")

dev.off()

# Make figure (Figure 2D-E) with both Boxplots. For Fig 2D, HOF1 is removed because the KO strain has almost no growth.The figure was also edited in illustrator.

f2d <- diff_with_WT %>% filter(!is.na(diff),gene_tagged !="HOF1")%>%
  ggplot(., aes(x=reorder(gene_tagged, diff,fun=median), y=diff, fill=type_diff)) +
  geom_boxplot(outlier.size=-1)+ 
  theme_bw()+
  ylim(-5,1)+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'))+
  scale_fill_manual(values=c("grey1","grey78"),name = "Comparison versus WT", labels = c("Deletion of SH3","Gene deletion"))+
  labs(y="Growth relative to WT", x="Gene")+
  coord_flip()+
  theme(legend.position = c(0.2, 0.8))


f2e <-ggplot(mSH3_ko_d_scores, aes(x=reorder(orf, desc(rank),fun=median), y=rel_to_wt, fill=type_construct))+
  geom_boxplot(outlier.size=-1)+
  theme_bw()+
  ylim(-5,1)+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'))+
  scale_fill_manual(values=c("grey1","grey78"),name = "Comparison versus WT", 
                    labels = c("Deletion of SH3","Gene deletion"))+
  labs(y="Growth relative to WT", x="Gene")+
  scale_x_discrete(labels = c("D|D|D","WT|WT|D","WT|D|WT","D|WT|WT","SLA1","D|D","WT|D","D|WT","BEM1","D|D","WT|D","D|WT","BZZ1") )+
  coord_flip()+
  theme(legend.position = "none")



pdf("BoxPlotStuffedKOVsWTALL_CRLCode_0320.pdf",7, 10)
plot_grid(f2d, f2e, labels = c('D', 'E'), 
          nrow=2,rel_heights =  c(20, 13),
          align = 'v')
dev.off()                              

# The analyzed datasets that are available with the article correspond to f2d and this one:
mSH3 = mSH3_ko_d_scores %>%  select(orf, gene_tagged, Condition, type_construct, rel_to_wt)
write.csv(mSH3,"GrowthMultipleSH3Analyzed_0420.csv")
# Some column names and strain names have been manually changed for the paper.