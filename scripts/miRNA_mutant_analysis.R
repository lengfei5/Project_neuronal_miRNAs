##########################################################################
##########################################################################
# Project: Chiara's neuronal project
# Script purpose: analysis the mutant data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 27 15:57:22 2019
##########################################################################
##########################################################################
version.Data = 'miRNAs_neurons_v1';
version.analysis = paste0(version.Data, "_2019_05_27")

miRNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/miRNAseq_functions.R"

### Directories to save results
design.file = "../exp_design/Neuron_project_design_all_v1.xlsx"
dataDir = "../data"
resDir = paste0("../results/miRNAs_neurons_mutants")
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##################################################
##################################################
## Section: Import Sample information and table of read counts
##################################################
##################################################
library("openxlsx")
design = read.xlsx(design.file, sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, skipEmptyCols = TRUE)

Processing.design.matrix = TRUE
if(Processing.design.matrix){
  xx = design;
  design = data.frame(xx$sample.ID, as.character(xx$genotype), xx$`Tissue/.Cell-type`, xx$promoter, xx$sampleInfo, xx$sampleInfo, 
                      stringsAsFactors = FALSE)
  colnames(design) = c('SampleID', 'genotype', 'tissue.cell', 'promoter', 'treatment', 'sampleInfo')
  
  kk = grep("No Treatment|No treatment", design$treatment)
  design$treatment[kk] = "untreated"
  kk = grep("Treatment", design$treatment)
  design$treatment[kk] = "treated"
  jj = grep("henn-1", design$genotype)
  design$genotype[jj] = "henn1.mutant"
  #design$genotype[grep("WT", design$genotype)] = "WT"
  
  design$promoter = sapply(design$promoter, function(x) unlist(strsplit(as.character(x), ":"))[1])
  design$promoter = sapply(design$promoter, function(x) gsub('_', '.', x))
  design$tissue.cell = sapply(design$tissue.cell, function(x) gsub(' ', '.', x))
  
  #design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole.body_no_promoter"
}

jj = c(which(design$genotype == "N2"), which(design$promoter == "osm-5"))
design = design[jj, ]


## make the data table
xlist <- list.files(path=paste0(dataDir), pattern = "*.txt", full.names = TRUE) ## list of data set to merge
if(length(xlist)>1){
  all = NULL
  ggs = NULL
  #counts = list()
  for(n in 1:length(xlist)){
    ff = read.delim(xlist[n], sep='\t', header = TRUE);
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
  all = read.delim(paste0(dataDir, "cel_countTable.txt"), sep = "\t", header = TRUE)
}

# processing count table for data before new nf-pipeline 
source(miRNAfunctions)
jj = c(which(design$genotype == "N2"))
all = process.countTable.v1(all=all, design = design[jj, ]); # process.countTable.v1 is the function version used before

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
aa = read.delim(paste0(dataDir, "/R7794_nf_all/countTable.txt"), sep = "\t", header = TRUE)

# spike-ins
spikes = read.delim(paste0(dataDir, "/R7794_nf_all/spikeInTable.txt"), sep="\t", header = TRUE, row.names = 1)

# start to process table and merge them into one
kk = grep("piRNA_", aa$Name)
if(length(kk)>0){
  piRNAs = aa[kk, ]
  aa = aa[-kk,]
}

source(miRNAfunctions)

jj2 = which(design$promoter == "osm-5")
all2 = process.countTable(all=aa, design = design[jj2, c(1:5)], select.counts = "Total.count")

xx = as.matrix(all[, -1])
xx[which(is.na(xx))] = 0
stat.miRNAs = floor(apply(xx, 2, sum))
piRNAs = process.countTable(all= piRNAs, design = design[jj2, c(1:5)], select.counts = "GM.count")

piRNAs = as.matrix(piRNAs[, -1])
piRNAs[which(is.na(piRNAs))] = 0
stat.piRNAs = floor(apply(piRNAs, 2, sum))
#all = rbind(c('total.piRNAs', stat.piRNAs), all)

spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE)
spikes = process.countTable(all=spikes, design =  design[jj2, c(1:5)], select.counts = "Total.spikeIn")

total.spikes = floor(apply(as.matrix(spikes[, -1]), 2, sum))

##########################################
# quick comparison between old pipeline and new pipeline from nf-all 
##########################################
Compare.old.and.new.nf_pipeline = FALSE
if(Compare.old.and.new.nf_pipeline){
  mm = match(all$gene, all2$gene)
  for(id2check in c("52632", "52633", "52634", "52635"))
  {
    yy1 = all[, grep(id2check, colnames(all))]
    yy2 = all2[mm, grep(id2check, colnames(all2))]
    plot(yy1, yy2, log='xy'); abline(0, 1, lwd=2.0, col='red')
    #abline(v=23)
    cat(" id : ", id2check, "--", all$gene[which(abs(yy1-yy2)>1)], "--diff --", (yy1-yy2)[which(abs(yy1-yy2)>1)] ,"..\n")
  }
  
  # check the piRNA counts
  load( file = paste0("../results/tables_for_decomvolution/Rdata/", 
                      'neuronalClasses_samples_countTables_piRAN_siRNA_stats_', 
                      "miRNAs_neurons_v1_2018_03_07", '.Rdata'))
  
  plot(stat.piRNAs[c(1:4)], stats$piRNA[grep("52632|52633|52634|52635", rownames(stats))], log='xy');abline(0, 1, lwd=2.0, col='red')
  
}

##########################################
#  combine the old data set and new data from nf-all
##########################################
all = rbind(spikes, all);

design.matrix = data.frame(design, stat.miRNAs, stat.piRNAs, total.spikes)
