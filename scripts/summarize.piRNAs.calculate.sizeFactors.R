##########################################################################
##########################################################################
# Project:
# Script purpose: to count individual piRNA for the normalization
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct  2 14:02:27 2018
##########################################################################
##########################################################################
analysisAllpiRNA <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", 
                             file.suffixSecondary="bam.seqCnt.txt.gz", 
                             includePrimary="NONE", 
                             includeType = c("piRNA"), 
                             tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE, countColumn = "count")
{
  source("functions.miRNASummarize.noFixation.R")
  dataDir = "../data/normalized_piRNAs/counts_smRNA"
  
  library(data.table)
  
  file.suffixPrimary = "bam.seqCnt.txt"; file.suffixSecondary="bam.seqCnt.txt.gz";
  includePrimary="NONE"; includeType = c("piRNA"); annotWithAS=TRUE; 
  normalizeCnt = FALSE; tailFilter = 0.12;
  excludeId = NULL; excludeType = NULL;
  lenDis = 18:30; aboveIsoCnt = 0; processes = 1; spliceSite = FALSE;
  last8mer = FALSE;  minLength=-1; countColumn="count";
  
  seqCnt.all <- list.files(dataDir, paste0(file.suffixPrimary, "$"), full.names = TRUE)
  
  for(n in 1:length(seqCnt.all))
  {
    # n = 1
    ff = fread(seqCnt.all[n])
    ff = ff[which(ff$type=="piRNA;" & ff$tailLen == 0), ]
    #ff = read.delim(seqCnt.all[n], sep = "\t")
  }
  #miRNA.list <- lapply(seqCnt.all, function(seqCnt.file)
  #{
  #  show(seqCnt.file)
  # seqCnt <- getSeqCnt(seqCnt.file=seqCnt.file, file.suffixPrimary=file.suffixPrimary, 
  #                      file.suffixSecondary=file.suffixSecondary, normalizeCnt=normalizeCnt, 
  #                      tailFilter=tailFilter, includePrimary=includePrimary, includeType=includeType, 
  #                      excludeId=excludeId, excludeType=excludeType, spliceSite=spliceSite, annotWithAS=annotWithAS, minLength=minLength, countColumn=countColumn)
  #  return(analyseSeqCnt(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer))
  #})
  
  
  #test.script = FALSE
  #if(test.script){
  #  }
  
  #$analysisAll(file.suffixPrimary = file.suffixPrimary, 
  #            output.middleName = "piRNA", 
  #            file.suffixSecondary = file.suffixSecondary, 
  #            tailFilter = tailFilter, includePrimary = includePrimary, 
  #            includeType = includeType, 
  #            aboveIsoCnt = aboveIsoCnt, 
  #            processes = processes, 
  #            annotWithAS=annotWithAS,
  #            countColumn = countColumn)
  
}
