dada2_calc_trims <- function(trim_low,trim_high,...){
  
  for(trim in trim_low:trim_high){
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(trim,trim),
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
    
    if(trim == trim_low){
      track_all <- track
    } else 
      track_all <- track_all %>% full_join(track)
  }
  return(track_all)
}