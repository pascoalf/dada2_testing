dada2_calc_trims_groups <- function(trims_F,trims_R,...){
require(dada2)
require(vegan)
  for(trim in 1:length(trims_F)){
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(trims_F[trim],trims_R[trim]),
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
    
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    seqtab <- makeSequenceTable(mergers)
    
    #this is the final ASV table
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    
    ##track reads
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
    
    rownames(track) <- sample.names
    track <- as.data.frame(track)
    
    #add number of ASVs into track table
    
    track$species_richness <- specnumber(seqtab.nochim)
    track$samples <- row.names(track)
    track$trimming <- trim
    row.names(track) <- NULL
    track <- as.data.frame(track)
    
    if(trim == 1){
      track_all <- track
    } else 
      track_all <- track_all %>% full_join(track)
  }
  return(track_all)
}