##########################################################################
##########################################################################
# Project:
# Script purpose: to count individual piRNA for the normalization
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct  2 14:02:27 2018
##########################################################################
##########################################################################
count.piRNA.with.design.matrix = function(design.matrix, recount = FALSE)
{
  if(recount){
    library(data.table)
    dataDir = "../data/normalized_piRNAs/counts_smRNA"
    seqCnt.all <- list.files(dataDir, paste0(file.suffixPrimary, "$"), full.names = TRUE)
    
    kk =  unique(c(which(design$tissue.cell=="Glial.cells" | design$tissue.cell == "CEPsh" | 
                           design$genotype == "henn1.mutant"), grep("L3", design$tissue.cell)))
    
    design.matrix = design[-kk, ]
    
    
    test.if.samples.there = FALSE
    all = NULL
    ggs = NULL
    ids = NULL
    #counts = list()
    for(n in 1:nrow(design.matrix))
      #for(n in 1:5)
    {
      cat("--", n, "--\n")
      # n = 1
      #start_time <- Sys.time()
      #sleep_for_a_minute()
      #end_time <- Sys.time()
      jj = grep(design.matrix$SampleID[n], seqCnt.all)
      
      if(length(jj)==1) {
        ids = c(ids, design.matrix$SampleID[n])
        if(!test.if.samples.there){
          ff = fread(seqCnt.all[jj], nThread = 6)
          ff = ff[which(ff$type=="piRNA;" & ff$tailLen == 0), ]
          ggs = unique(ff$gname)
          
          counts = sapply(ggs, function(x) {sum(ff$count[which(ff$gname==x)])})
          
          ff = data.frame(ggs, counts, stringsAsFactors = FALSE)
          colnames(ff)[1] = "gene";
          
          if(n == 1){
            all = data.frame(ff, stringsAsFactors = FALSE);
            colnames(all)[1] = "gene"
          }else{
            ggs = ff[,1];
            gg.union = unique(union(all$gene, ggs))
            all = data.frame(gg.union,
                             all[match(gg.union, all$gene), -1],
                             ff[match(gg.union, ggs), -1], 
                             stringsAsFactors = FALSE)
            colnames(all)[1] = "gene"
          }
        }
        
      }else{
        
        cat(design.matrix$SampleID[n],  "-- NOT Found --\n")
      }
      
      #Sys.time() - start_time 
    }
    colnames(all)[-1] = ids;
    
    piRNA.counts = all;
    save(piRNA.counts, file = "../results/tables_for_decomvolution/Rdata/rawCounts_table_for_piRNAs.Rdata")
  }else{
    load(file = "../results/tables_for_decomvolution/Rdata/rawCounts_table_for_piRNAs.Rdata")
  }
  
  
}


##########################################
# try to use Thomas' function to count piRNAs but give up afterwards
##########################################
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
