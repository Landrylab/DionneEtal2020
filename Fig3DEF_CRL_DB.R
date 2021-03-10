library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot) #to use with plot_grid()
library(growthcurver)
library(ggrepel)
library(forcats)
library(seqinr)
library(ggExtra)
library(cowplot)
library(grid)
library(gridExtra)
rm(list=ls())


# Obtain the DMS data
dms <- read.table("SH3paperData - Fig4.DMSData.csv", sep=",", header=T)
dms_codon <- read.table("SH3paperData - Fig4.DMSData.csv", sep=",", header=T)

# Keep only HUA2_MTX2 and LSB3_MTX2
conditions<- c("HUA2_MTX2","LSB3_MTX2")

# Keep only the data for ABP1 domain 
dms %<>% filter(condition %in% conditions, sample=="ABP1", mutation.format=="amino acid")
dms_codon %<>% filter(condition %in% conditions, sample=="ABP1", mutation.format=="codon")

# Separate columns by amino acid and positon and rename some columns
dms %<>% separate(position.reference,c("pos", "wtaa"), sep=" ") %>%
        mutate(pos=as.integer(pos)+1) %>%
        rename(dms_score=ratio..log2.scale..mean) %>%
        select(sample, condition, mutation.format, pos, wtaa, mutated, dms_score,mutation.type)

# Do the same for data at the codon level
dms_codon %<>% separate(position.reference,c("pos", "wtcod"), sep=" ") %>%
  mutate(pos=as.integer(pos)+1) %>%
  rename(dms_score=ratio..log2.scale..mean) %>%
  select(sample, condition, mutation.format, pos, wtcod, mutated, dms_score,mutation.type)

# Find synonymous codons

library(Biostrings)
GENETIC_CODE

codon_mutations <- dms_codon[,c(5,6)]
status_vec <- NULL

# Designate codons as either synonymous, non-synonymous, or stop

for (i in 1:nrow(codon_mutations)) {
  
  print(i)
  
  codon_ref <- unname(codon_mutations[i,1])
  codon_mut <- unname(codon_mutations[i,2])
  
  aa_ref <- unname(GENETIC_CODE[names(GENETIC_CODE) == codon_ref])
  aa_mut <- unname(GENETIC_CODE[names(GENETIC_CODE) == codon_mut])
  
  if (aa_mut == '*') {status_vec <- c(status_vec, 'stop'); next}
  if (aa_ref == aa_mut) {status_vec <- c(status_vec, 'synonymous')}
  if (aa_ref != aa_mut) {status_vec <- c(status_vec, 'non-synonymous')}
  
}

dms_codon[,8] <- status_vec

# Have the dms scores of two PPIs side by side
datw <- dms %>% pivot_wider(names_from=condition, values_from = dms_score)

# Scale the data to have the synonymous mutations around 0
median_syn_codon <- dms_codon %>% group_by(condition,mutation.type) %>%
  summarize(median_score = median(dms_score, na.rm=T))

datw$hua2sc <- datw$HUA2_MTX2 - median_syn_codon$median_score[3]
datw$lsb3sc <- datw$LSB3_MTX2 - median_syn_codon$median_score[6]

# Add minimum distance of each SH3 domain residue to the peptide
d <- read_delim("2rpn_interchain_min_dist.txt", delim = " ")
d <- d[,c(1,3)]
colnames(d) <- c('pos','dis_pep')

#d %<>% rename(pos = ABP1, dis_pep = min_dist) %>% select(pos,dis_pep)
datw <- left_join(datw, d, by=c("pos"="pos"))

# Retrieve RSA values
rsa <- read_delim("abp1_2rpn_RSA.txt",col_names = FALSE, delim=" ")
names(rsa) <- c("pos", "rsa")
datw <- left_join(datw, rsa, by=c("pos"="pos"))


# Any mutation is considered deleterious if it falls below the 1st percentile
# of synonymous codon mutations. Calculate separately for HUA2 and LSB3

syn_hua2 <- dms_codon[dms_codon[,2] == 'HUA2_MTX2'&dms_codon[,8] == 'synonymous',7]
syn_hua2_sc <- syn_hua2 - median_syn_codon$median_score[3]

syn_lsb3 <- dms_codon[dms_codon[,2] == 'LSB3_MTX2'&dms_codon[,8] == 'synonymous',7]
syn_lsb3_sc <- syn_lsb3 - median_syn_codon$median_score[6]

# 1st percentiles  
q_syn_hua2 <- quantile(syn_hua2_sc, probs=seq(0,1,0.01),na.rm=T)[2]
q_syn_lsb3 <- quantile(syn_lsb3_sc, probs=seq(0,1,0.01),na.rm=T)[2]

# Assign each mutation to one of 4 possible quadrants

datw %<>% mutate(quad = ifelse(hua2sc < q_syn_hua2 & lsb3sc >= q_syn_lsb3,"Dest. Hua2", 
                        ifelse(hua2sc >= q_syn_hua2 & lsb3sc >= q_syn_lsb3,"Neutral",
                        ifelse(hua2sc < q_syn_hua2 & lsb3sc < q_syn_lsb3,"Dest. both",
                        ifelse(hua2sc >= q_syn_hua2 & lsb3sc < q_syn_lsb3,"Dest. Lsb3","NA")))))       

# For each quadrant, calculate the total number of mutations and number of SH3 positions (/58)

spec_to_lsb3 <- length(datw$pos[datw$quad=="Dest. Lsb3" & datw$mutation.type=="non-synonymous"])
nres_lsb3 <- length(unique(datw$pos[datw$quad=="Dest. Lsb3"& datw$mutation.type=="non-synonymous"]))

spec_to_hua2 <- length(datw$pos[datw$quad=="Dest. Hua2"& datw$mutation.type=="non-synonymous"])
nres_hua2 <- length(unique(datw$pos[datw$quad=="Dest. Hua2"& datw$mutation.type=="non-synonymous"]))

affect_both <- length(datw$pos[datw$quad=="Dest. both"& datw$mutation.type=="non-synonymous"])
nres_affect_both  <- length(unique(datw$pos[datw$quad=="Dest. both"& datw$mutation.type=="non-synonymous"]))

neutral_both <- length(datw$pos[datw$quad=="Neutral"& datw$mutation.type=="non-synonymous"])
nres_neutral <- length(unique(datw$pos[datw$quad=="Neutral"& datw$mutation.type=="non-synonymous"]))

# Calculate the correlation for the two DMS scores

test_corr <- cor.test(datw$lsb3sc, datw$hua2sc, method="kendall")

test_corr$p.value
test_corr$estimate

label_plot = paste("r = ", round(test_corr$estimate,2), ", p-value = 2.02e-205 ", sep="")

# Plot the correlation between HUA2 DMS scores and LSB3 DMS scores

plot_correlation <- 
  datw %>% rename(distance.peptide = dis_pep) %>%
  ggplot(.)+
  geom_point(aes(x=hua2sc, y=lsb3sc, shape=mutation.type, color=distance.peptide), alpha=0.9, size=3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'),
        legend.title = element_text(size = 11.5),
        legend.text = element_text(size = 10.5))+
  #scale_color_gradient(low=colors()[258], high=colors()[86])+
  xlab('DMS Score Abp1-Hua2')+
  ylab('DMS Score Abp1-Lsb3')+
  labs(shape='Mutation type')+
  labs(color='Peptide distance')+
  theme(legend.position="top", legend.box = "horizontal")+
  #guides(color=guide_legend(title.position = "left"))+
  guides(color=guide_colorbar(order=1))+
  theme(legend.key = element_rect(fill = "white")) +
  xlim(-3.2,2.5)+
  ylim(-3.2,2.5)+
  geom_abline(intercept = 0, slope = 1, color="grey50", linetype=3)+
  annotate(geom="text", x=1, y=2.5, label=label_plot,
           size=5,col='black')+
  geom_vline(xintercept = q_syn_hua2, color="grey50", linetype=3)+
  geom_hline(yintercept = q_syn_lsb3, color="grey50", linetype=3)+
  annotate(geom="text", x=-2.25, y=2.5, 
           label=paste("Destabilizing for Abp1-Hua2", "\nn=", spec_to_hua2, 
                       ", ", nres_hua2, " positions", 
                       sep=""),
           size=4.5)+
  annotate(geom="text", x=0.6, y=-3.2, 
           label=paste("Destabilizing for Abp1-Lsb3", "\nn=", spec_to_lsb3, 
                       ", ", nres_lsb3, " positions",
                       sep=""),
           size=4.5)+
  annotate(geom="text", x=-2,, y=-3.2, 
           label=paste("Destabilizing for both", "\nn=", affect_both, 
                       ", ", nres_affect_both, " positions",
                       sep=""),
           size=4.5)


# Boxplot of distances for the 4 quadrants

boxplot_distance <- 
  datw %>% rename(distance.peptide = dis_pep) %>% 
  filter(mutation.type=="non-synonymous", !is.na(quad)) %>%
  ggplot(.)+
  geom_boxplot(aes(x=quad, y=distance.peptide), fill="grey",lwd=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'),
        legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 12))+
  xlab('')+
  ylab('Distance from peptide (Å)')+ 
  coord_flip()

boxplot_distance <- boxplot_distance + theme(plot.margin = unit(c(1.2,0.7,0,0), "cm")) #t, r, b, l

# Boxplot of RSAs for the 4 quadrants

boxplot_rsa <- 
  datw %>% rename(distance.peptide = dis_pep) %>% 
  filter(mutation.type=="non-synonymous", !is.na(quad)) %>%
  ggplot(.)+
  geom_boxplot(aes(x=quad, y=rsa), fill="grey",lwd=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13.5,face="bold"),
        plot.title=element_text(size=14,face='bold'),
        legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 12))+
  xlab('')+
  ylab('RSA')+
  coord_flip()

boxplot_rsa <- boxplot_rsa + theme(plot.margin = unit(c(0.9,0.7,0.3,0), "cm")) #t, r, b, l

# Assemble the 3 plots now using Cowplot

lower_plot <- plot_grid(boxplot_distance, boxplot_rsa,nrow=2,align = 'v')
y.grob <- textGrob('Type of non-synonymous mutation', 
                   gp=gpar(fontface="bold", fontsize=13.5),rot=90)

lower_plot <- grid.arrange(arrangeGrob(lower_plot, left = y.grob))
dms_figure <- plot_grid(plot_correlation, lower_plot,nrow=1,
                        rel_widths = c(1.7,1))

dms_figure

# Plot in RStudio with PDF, A4, landscape
