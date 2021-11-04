##dada2 extras
#data from tutorial
library(dada2); packageVersion("dada2")
library(dplyr)
library(vegan)
library(ggplot2)
#
path <- "./samples"
#
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#

#inspect read quality

pdf("./figures/quality profiles of forward reads.pdf")
plotQualityProfile(fnFs[1:length(fnRs)])
dev.off()
pdf("./figures/quality profiles of reverse reads.pdf")
plotQualityProfile(fnFs[1:length(fnRs)])
dev.off()

#
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#this function runs dada2 from the trimming step to the final ASV table
#no taxonomic assignments
source("dada2_calc_trims.r")

##test maxEE
#this chunk takes 43 minutes in linux computer
start_1 <- Sys.time()
track_test_maxE_2.2 <- dada2_calc_trims(trim_low = 200, trim_high = 201) %>%  #maxEE is the custom
                          mutate(maxEE = "2.2", rm.phix=T)
end_1 <- Sys.time()
start_2 <- Sys.time()
track_test_maxE_2.3 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,3)) %>% mutate(maxEE = "2.3", rm.phix=T)
end_2 <- Sys.time()
start_3 <- Sys.time()
track_test_maxE_2.4 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,4)) %>% mutate(maxEE = "2.4", rm.phix=T)
end_3 <- Sys.time()
start_4 <- Sys.time()
track_test_maxE_2.5 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,5)) %>% mutate(maxEE = "2.5", rm.phix=T)
end_4 <- Sys.time()
start_5 <- Sys.time()
track_test_maxE_3.2 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(3,2)) %>% mutate(maxEE = "3.2", rm.phix=T)
end_5 <- Sys.time()
start_6 <- Sys.time()
track_test_maxE_4.2 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(4,2)) %>% mutate(maxEE = "4.2", rm.phix=T)
end_6 <- Sys.time()
start_7 <- Sys.time()
track_test_maxE_5.2 <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(5,2)) %>% mutate(maxEE = "5.2", rm.phix=T)
end_7 <- Sys.time()

##test phy

start_8 <- Sys.time()
track_test_maxE_2.2_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,rm.phix=F) %>%  mutate(maxEE = "2.2", rm.phix=F)
end_8 <- Sys.time()
start_9 <- Sys.time()
track_test_maxE_2.3_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,3),rm.phix=F) %>% mutate(maxEE = "2.3", rm.phix=F)
end_9 <- Sys.time()
start_10 <- Sys.time()
track_test_maxE_2.4_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,4),rm.phix=F) %>% mutate(maxEE = "2.4", rm.phix=F)
end_10 <- Sys.time()
start_11 <- Sys.time()
track_test_maxE_2.5_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(2,5),rm.phix=F) %>% mutate(maxEE = "2.5", rm.phix=F)
end_12 <- Sys.time()
start_13 <- Sys.time()
track_test_maxE_3.2_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(3,2),rm.phix=F) %>% mutate(maxEE = "3.2", rm.phix=F)
end_13 <- Sys.time()
start_14 <- Sys.time()
track_test_maxE_4.2_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(4,2),rm.phix=F) %>% mutate(maxEE = "4.2", rm.phix=F)
end_14 <- Sys.time()
start_15 <- Sys.time()
track_test_maxE_5.2_nophy <- dada2_calc_trims(trim_low = 200, trim_high = 201,maxEE = c(5,2),rm.phix=F) %>% mutate(maxEE = "5.2", rm.phix=F)
end_15 <- Sys.time()

###
track_test_maxE_2.3_nophy<-track_test_maxE_2.3_nophy %>% mutate(rm.phix = F)
track_test_maxE_2.4_nophy<-track_test_maxE_2.4_nophy %>% mutate(rm.phix = F)
track_test_maxE_2.5_nophy<-track_test_maxE_2.5_nophy %>% mutate(rm.phix = F)
track_test_maxE_3.2_nophy<-track_test_maxE_3.2_nophy %>% mutate(rm.phix = F)
track_test_maxE_4.2_nophy<-track_test_maxE_4.2_nophy %>% mutate(rm.phix = F)
track_test_maxE_5.2_nophy<-track_test_maxE_5.2_nophy %>% mutate(rm.phix = F)
###

##join tables for figures
track_combined <- rbind(track_test_maxE_2.2,
                        track_test_maxE_2.3,
                        track_test_maxE_2.4,
                        track_test_maxE_2.5,
                        track_test_maxE_3.2,
                        track_test_maxE_4.2,
                        track_test_maxE_5.2,
                        track_test_maxE_2.2_nophy,
                        track_test_maxE_2.3_nophy,
                        track_test_maxE_2.4_nophy,
                        track_test_maxE_2.5_nophy,
                        track_test_maxE_3.2_nophy,
                        track_test_maxE_4.2_nophy,
                        track_test_maxE_5.2_nophy)
#
pdf("amostra de resultados - dada2 testing.pdf")
track_combined %>% 
  ggplot(aes(x=as.factor(trimming),y=nonchim,col=rm.phix))+
  geom_boxplot()+facet_wrap(~ maxEE)+
  labs(title = "Boxplot for number of sequences, for several trimmings, maxEE and rm.phix")
#
track_combined %>% ggplot(aes(x=as.factor(trimming),y=(input-nonchim)*100/input,col=rm.phix))+
  geom_boxplot()+facet_wrap(~ maxEE)+expand_limits(y=c(0,100))+
  labs(title = "Boxplot for % of lost sequences, for several trimmings, maxEE and rm.phix")
dev.off()

### check the time for one full trimming
## start 1 to end 1 was the trimming of 200:201
## calculate the time for 1 more nt
## i.e. 200:202

start_202 <- Sys.time()
testing_202 <- dada2_calc_trims(trim_low = 200, trim_high = 202) %>% mutate(maxEE = "2.2",rm.phix = T)
end_202 <- Sys.time()

##Check if is possible to do one full trimming

start_full <- Sys.time()
testing_full <- dada2_calc_trims(trim_low = 190, trim_high = 250) %>% 
  mutate(maxEE = "2.2",rm.phix = T)
end_full <- Sys.time()

