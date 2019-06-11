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
QCfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_Quality_Controls.R"

### Directories to save results
design.file = "../exp_design/Neuron_project_design_all_v1.xlsx"
dataDir = "../data"
resDir = paste0("../results/miRNAs_neurons_mutants")
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

Compare.old.and.new.nf_pipeline = FALSE

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

all2 = rbind(spikes, all2);
##########################################
# quick comparison between old pipeline and new pipeline from nf-all 
##########################################
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
xx = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(all, all2))
#design.matrix = data.frame(design, stat.miRNAs, stat.piRNAs, total.spikes)
all = xx
save(design, all, stat.piRNAs, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

########################################################
########################################################
# Section : merge technical replicates and filter genes
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
source(miRNAfunctions)

Filter.lowly.expressed.using.predefined.miRNA.list = TRUE
Filter.lowly.expressed.using.cpm.threshold = FALSE
Merge.techinical.replicates.N2 = TRUE

tcs = unique(design$tissue.cell)
tcs = setdiff(tcs, c("whole.body", "whole.body.L3"))
length(tcs)

if(Merge.techinical.replicates.N2){
  rep.technical = list(c("57751", "57753"), c("57752", "57754"))
  res.merged = Merge.techinical.replicates.using.design.countTable(design, all, id.list=rep.technical);
  all = res.merged$countTable
  design = res.merged$design
}

# filter lowly expressed miRNA with list of predefined miRNAs that were identified using all untreated samples 
if(Filter.lowly.expressed.using.predefined.miRNA.list){
  list.expressed = read.csv(paste0(dataDir, "/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv"), 
                            header = TRUE, as.is = c(1, 2))
  expressed.miRNAs = find.mature.ones.for.prefixed.expressed.miRNAs(list.expressed)
}

if(Filter.lowly.expressed.using.cpm.threshold){ # filter the non-expressed miRNAs using cpm threshold
  
  kk = intersect(grep("L3", design$tissue.cell), which(design$treatment=='untreated'))
  rownames(all) = all$gene 
  ## the rownames should be the miRNA names
  expressed.miRNAs = find.expressed.mature.miRNA.using.cpm.threshold(all[, (kk+1)], cpm.threshold=10)
  
}

##########################################
# check the quality using DESeq2 
##########################################
read.count = all[, -1]
kk = c(1:nrow(design))
design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
raw = as.matrix(read.count[,kk])
raw[which(is.na(raw))] = 0
### start enrichment analysis 
raw = floor(raw)
rownames(raw) = all$gene

index.qc = c(1, 3, 6)

#norms = as.numeric(res.spike.in$norms4DESeq2)
#norms = as.numeric(piRNAs)/median(as.numeric(piRNAs))

source(QCfunctions)

pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, ".pdf")

pdf(pdfname, width = 18, height = 10)
Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix[, index.qc])
dev.off()

##########################################
# Normalize the all data with DESeq2  
##########################################
dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ genotype + treatment)
#dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
dds1 <- dds[ rowSums(counts(dds)) > 10, ]
dds1 <- estimateSizeFactors(dds1)
sizeFactors(dds) = sizeFactors(dds1)
cpm = fpm(dds, robust = TRUE)
colnames(cpm) = paste0(colnames(cpm), ".DESeq2norm")
xx = log2(cpm + 0.1)

pdfname = paste0(resDir, "/pairwise_comparisons_untreated_DESeq2_", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 16)
pairs(xx[, grep('_untreated', colnames(xx))], lower.panel=NULL, upper.panel=panel.fitting, cex = 0.5)
dev.off()

##########################################
# spike-in normalization
##########################################
index.spikeIn = grep("spikeIn", rownames(raw))[c(1:8)]
spike.concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100 ## the concentration is amol/mug-of-total-RNA

## calculate scaling factor using spike-ins
source(miRNAfunctions)

jj = grep("WT_tax4", design$genotype)
pdfname = paste0(resDir, "/Spike_in_signals_normalized_", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 10)
par(mfrow=c(2,2))

res.spike.in = calculate.scaling.factors.using.spikeIns(raw[,jj], concentrations = spike.concentrations, index.spikeIn = index.spikeIn, read.threshold = 5)

dev.off()

#cpm = res.spike.in$cpm;
res = res.spike.in$normalization.spikeIn
#colnames(cpm) = paste0(colnames(cpm), ".cpm")
colnames(res) = paste0(colnames(res), ".amol.per.mugRNA.normBySpikeIns")

pdfname = paste0(resDir, "/pairwise_comparisons_spikeIns_normalization_untreated_", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 12)
xx = log2(res + 0.01)
pairs(xx[, grep('_untreated', colnames(xx))], lower.panel=NULL, upper.panel=panel.fitting, cex = 0.7)
dev.off()

##########################################
# piRNA normalization
##########################################
# check the piRNA counts
load( file = paste0("../results/tables_for_decomvolution/Rdata/", 
                    'neuronalClasses_samples_countTables_piRAN_siRNA_stats_', 
                    "miRNAs_neurons_v1_2018_03_07", '.Rdata'))

piRNAs = c()
for(id in design$SampleID[c(1:6)]){
  piRNAs = c(piRNAs, stats$piRNA[grep(id, rownames(stats))])
}

piRNAs = c(piRNAs, stat.piRNAs)
#piRNAs = design.matrix$stat.piRNAs
sizefactors = as.numeric(piRNAs)
cpm.piRNA = raw
for(n in 1:ncol(cpm.piRNA))
{
  cpm.piRNA[,n] = raw[,n]/sizefactors[n]*10^6
}
colnames(cpm.piRNA) = paste0(colnames(cpm.piRNA), "normBy.piRNA")

pdfname = paste0(resDir, "/pairwise_comparisons_piRNAs_normalization_untreated_", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 12)
xx = log2(cpm.piRNA + 0.01)
pairs(xx[, grep('_untreated', colnames(xx))], lower.panel=NULL, upper.panel=panel.fitting, cex = 0.4)
dev.off()

if(Save.Tables){
  yy = data.frame(cpm, res, cpm.piRNA, stringsAsFactors = FALSE)
  write.csv(yy, file = paste0(tabDir, "Normalized_Table_DESeq2_spikeIn_piRNAs_normalized_for", version.analysis, ".csv"), 
            row.names = TRUE)
}

#plot(xx[, c(15, 16)]); abline(0, 1, lwd=2.0, col='red')

