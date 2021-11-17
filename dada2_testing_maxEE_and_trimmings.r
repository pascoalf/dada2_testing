##dada2 extras
#data from tutorial
library(dada2); packageVersion("dada2")
library(dplyr)
library(vegan)
library(ggplot2)
#please add dada2_calc_trims.r function
source("dada2_calc_trims.r")
#
path <- "./samples"
#
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

##test maxEE
start_1 <- Sys.time()
track_test_maxE_2.2 <- dada2_calc_trims(trim_low = 190, trim_high = 250) %>%  #maxEE is the custom
  mutate(maxEE = "2.2", rm.phix=T)
write.table(track_test_maxE_2.2,"track_test_maxE_2.2")
end_1 <- Sys.time()
start_2 <- Sys.time()
track_test_maxE_2.3 <- dada2_calc_trims(trim_low = 190, trim_high = 250,maxEE = c(2,3)) %>% mutate(maxEE = "2.3", rm.phix=T)
write.table(track_test_maxE_2.3,"track_test_maxE_2.3")
end_2 <- Sys.time()
start_3 <- Sys.time()
track_test_maxE_2.4 <- dada2_calc_trims(trim_low = 190, trim_high = 250,maxEE = c(2,4)) %>% mutate(maxEE = "2.4", rm.phix=T)
write.table(track_test_maxE_2.4,"track_test_maxE_2.4")
end_3 <- Sys.time()
##join tables for figures
track_combined <- rbind(track_test_maxE_2.2,
                        track_test_maxE_2.3,
                        track_test_maxE_2.4)
#
pdf("Testing some variables of dada2")
#
track_combined %>% 
  ggplot(aes(x=as.factor(trimming),y=nonchim))+
  geom_boxplot()+facet_wrap(~ maxEE)+
  xlab("Trimmings (190 to 250)")+
  ylab("Number of sequences")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  labs(title = "Boxplot for number of sequences, for several trimmings, maxEE")
#
track_combined %>% 
  ggplot(aes(x=as.factor(trimming),y=species_richness))+
  geom_boxplot()+facet_wrap(~ maxEE)+
  xlab("Trimmings (190 to 250)")+
  ylab("Species richness")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  labs(title = "Boxplot for number of sequences, for several trimmings, maxEE")
#
track_combined %>% ggplot(aes(x=as.factor(trimming),y=(input-nonchim)*100/input))+
  geom_boxplot()+facet_wrap(~ maxEE)+expand_limits(y=c(0,100))+ 
  xlab("Trimmings (190 to 250)")+
  ylab("Percentage of lost sequences")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  labs(title = "Boxplot for % of lost sequences, for several trimmings, maxEE")
#
dev.off()

####
source("dada2_calc_trims_groups.r")
f_trims <- c(240,220,200)
r_trims <- c(220,200,190)
start_groups_test <- Sys.time()
groups_test <- dada2_calc_trims_groups(trims_F = f_trims, trims_R = r_trims)
end_groups_test <- Sys.time()
groups_test <- groups_test %>% mutate(trimming_groups = paste(f_trims[trimming],r_trims[trimming]))

groups_test %>% ggplot(aes(x=as.factor(trimming_groups),y=nonchim))+geom_boxplot()+theme_classic()+
  labs(title="Total Sequences for each group of F and R trimmings")+ 
  xlab("Trimmings (F and R)") +
  geom_jitter(alpha=0.2)
  
groups_test %>% ggplot(aes(x=as.factor(trimming_groups),y=species_richness))+geom_boxplot()+theme_classic()+
  labs(title="Total Species richness for each group of F and R trimmings")+ 
  xlab("Trimmings (F and R)") +
  geom_jitter(alpha=0.2)
