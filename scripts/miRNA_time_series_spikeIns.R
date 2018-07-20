##################################################
## Project: Chiara's time serie miRNAs 
## Script purpose: control data quality, normalize data using spike-ins, enrichment analysis 
## the script was modified based on the spike-in normalization script for Philipp and Xue 
## The enrichment analysis and also time series analysis are also added in this script
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
library("openxlsx")
require('DESeq2')
source('miRNAseq_functions.R')

### data verision and analysis version   
version.Data = 'miRNAs_timeSeries_spikeIn_R5922_R6016';
version.analysis = paste0("_", version.Data, "_20180719")

### Directories to save results 
design.file = "../exp_design/Libaries_time_series_spikIns.xlsx"
dataDir = "../data/timeSeries_withSpikIns"

#data.file = paste0(dataDir, "countTable.txt")
#spikes.file = paste0(dataDir, "spikeIns_count_table.txt")

resDir = "../results/time_series"
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##################################################
##################################################
## Section: import design matrix and prepare count table
##################################################
##################################################
design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
#design = design[which(design$stage != "L1.early"),]

Processing.design.matrix = TRUE
if(Processing.design.matrix){
  xx = design;
  design = data.frame(xx$ID, xx$treatment, xx$stage, as.character(xx$background), xx$`Tissue/Cell-type`, xx$Sample, stringsAsFactors = FALSE)
  colnames(design) = c('SampleID',  'treatment', 'stage', 'genotype', 'tissue.cell', 'sampleInfo')
  
  kk = grep("No Treatment|No treatment", design$treatment)
  design$treatment[kk] = "untreated"
  kk = grep("Treatment", design$treatment)
  design$treatment[kk] = "treated"
  design = design[order(design$stage, design$treatment), ]
  #design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole.body_no_promoter"
}

## make the data table
xlist<-list.files(path=paste0(dataDir), pattern = "*.txt", full.names = TRUE) ## list of data set to merge
spikes.file = xlist[grep("spikeIns_", xlist)]
spikes.file = spikes.file[grep("_old", spikes.file, invert = TRUE)]
xlist = xlist[grep("cel_", xlist)]

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
  all = read.delim(xlist, sep = "\t", header = TRUE)
}
# processing count table
all = process.countTable(all=all, design = design);

#### Import Sample information and table of read counts
#design = read.xlsx(paste0(design.file), sheet = 1, colNames = TRUE)
#colnames(design)[c(1:2)] = c("SampleID", "genotype")
#design$genotype[grep("WT", design$genotype)] = "WT"
#design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole_animal_no_promoter"

#all = read.delim(data.file, sep = "\t", header = TRUE)
# processing count table
#all = process.countTable(all=all, design = design);

spikes = read.delim(spikes.file, sep="\t", header = TRUE, row.names = 1)
spikes = t(as.matrix(spikes))
spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE)
spikes = process.countTable(all=spikes, design = design, select.Total.count = FALSE)

all = rbind(spikes, all);
#spikes = process.countTable(spikes, design)

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

##################################################
##################################################
## Section: spike-in normalization and QC for cpm and spike-in normalization
##################################################
##################################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
read.count = all[, -1];
kk = c(1:nrow(design))

design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
raw = as.matrix(read.count[,kk])
raw[which(is.na(raw))] = 0
### start enrichment analysis 
raw = floor(raw)
rownames(raw) = all$gene

####################
## QC for cpm 
####################
QC.for.cpm = FALSE
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  index.qc = c(1, 3,4)
  
  source("RNAseq_Quality_Controls.R")
  pdfname = paste0(resDir, "/Data_qulity_assessment_early_L1_and_all", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix[, index.qc])
  dev.off()
}

####################
## spike-in normalization 
####################
sel.samples.with.spikeIns = which(design.matrix$stage != "L1.early")
design.matrix = data.frame(sample=colnames(read.count)[sel.samples.with.spikeIns], design[sel.samples.with.spikeIns, ])
countData = raw[, sel.samples.with.spikeIns]

dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ treatment + stage)
#dds$genotype <- relevel(dds$genotype, ref="WT");
#dds$treatment = relevel(dds$treatment, ref="untreated");
index.spikeIn = grep("spikeIn", rownames(dds))
concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100

## calculate scaling factor using spike-ins
source("miRNAseq_functions.R")
pdfname = paste0(resDir, "/Spike_in_signals_normalized_DESeq", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 10)

par(mfrow=c(2,2))
#norms = calculate.scaling.factors.using.spikeIns(dds, concentrations = concentrations, index.spikeIn = index.spikeIn, read.threshold = 5)
res.spike.in = calculate.scaling.factors.using.spikeIns(countData = countData, 
                                                        concentrations = concentrations, 
                                                        index.spikeIn = index.spikeIn, read.threshold = 5)

dev.off()

cpm = res.spike.in$cpm;
res = res.spike.in$normalization.spikeIn
colnames(cpm) = paste0(colnames(cpm), ".cpm")
colnames(res) = paste0(colnames(res), ".amol.per.mugRNA.normBySpikeIn")

norms = res.spike.in$norms4DESeq2;
sizeFactors(dds) = res.spike.in$norms4DESeq2

#res = fpm(dds, robust = TRUE)
#cpm =  fpm(dds, robust = FALSE)
#kk = grep("spikeIn", rownames(cpm))

ss = apply(countData, 2, sum)

plot(countData[,1]/ss[1]*10^6, cpm[,1], log='xy');abline(0, 1, lwd=2.0, col='red')
#plot(raw[,1]/ss[1]*10^6/norms[1], res[,1], log='xy');abline(0, 1, lwd=2.0, col='red')

#plot(raw[,1]/prod(ss)^(1/length(ss))*10^6/si, res[,1], log='xy');abline(0, 1, lwd=2.0, col='red')
res = data.frame(res, cpm)

write.table(res, file = paste0(tabDir, paste0("normalized_signals_scalling.factors_using_spikeIns", version.analysis, ".txt")), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

save(all, design, res, file=paste0(RdataDir, 'Design_Raw_readCounts_cpm_spikeInsNorm', version.analysis, '.Rdata'))

####################
## QCs for spike-in normalization 
####################
index.qc = c(1, 3,4)
source("RNAseq_Quality_Controls.R")
pdfname = paste0(resDir, "/Data_qulity_assessment_all_samples_withSpikeIns", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)
Check.RNAseq.Quality(read.count=countData, design.matrix = design.matrix[, index.qc], norms = norms)

dev.off()

####################
## Compare two versions of spike-in normalization 
####################
Compare.spikeIns.old.vs.new = FALSE
if(Compare.spikeIns.old.vs.new)
{
  version.analysis.old = "_miRNAs_timeSeries_spikeIn_R5922_R6016_20180503"
  load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis.old, '.Rdata'))
  old = all[c(1:8), -c(1:5)]
  
  load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
  new = all[c(1:8), -c(1:5)]
  
  pdfname = paste0(resDir, "/SpikeIns_comparison_new_vs_old", version.analysis, ".pdf")
  pdf(pdfname, width = 16, height = 16)
  
  par(mfrow=c(4,4))
  for(n in 1:ncol(new))
  {
    plot(old[,n]+0.5, new[, n], type='p', cex=2.0, pch=16, log='xy', col= 'darkblue', xlab='spikeIn read counts (old)', ylab='spikeIn read counts (new)')
    abline(0, median(new[,n]/old[,n]), lwd=2.0, col='red', untf = TRUE)
  }
  
  dev.off()
}

####################
## make a plot for enriched miRNAs
####################
list.enriched = read.csv("../results/miRNAs_neurons_v1_2018_03_07/tables/Pan.neurons/miRNA_Enrichment_Analysis_Pan.neurons_henn1.mutant_rab-3_Mature_miRNAs_neurons_v1_2018_03_07.csv", 
                         header = TRUE, row.names = 1)
list.enriched = list.enriched[order(-list.enriched$log2FoldChange), ]

pdfname = paste0(resDir, "/Expression_for_enriched_miRNAs_normalized_using_Spike_in_", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 8)

par(mfrow=c(1,1))
for(n in 1:nrow(list.enriched))
{
  #n = 1
  jj = grep(rownames(list.enriched)[n], rownames(res))
  #test = res[jj, ]
  
  #index = 1
  lty = 1
  if(length(jj)>1){
    index.untreated = grep("_untreated", colnames(res))
    ss.untreated = c()
    for(j in jj) ss.untreated = c(ss.untreated, mean(res[j, index.untreated]))
    jj = jj[which(ss.untreated==max(ss.untreated))][1] 
  }
  lims = range(res[jj,]) + 2^-4
  
  kk1 =  grep("_untreated", colnames(res))
  test = as.matrix(cbind(res[jj, kk1[which(c(1:length(kk1))%%2 ==1)]], res[jj, kk1[which(c(1:length(kk1))%%2 ==0)]]))
  matplot(c(1:4), (test+2^-4), col = 'darkblue', type = "p", lty = lty, ylim = lims, log='y', pch=1, main = rownames(res)[jj])
  points(c(1:4), apply(test, 1, mean)+2^-4, col='darkblue', type='l', lwd=2.0, lty= lty)  
  
  kk2 = grep("_treated", colnames(res))
  test = as.matrix(cbind(res[jj, kk2[which(c(1:length(kk1))%%2 ==1)]], res[jj, kk2[which(c(1:length(kk1))%%2 ==0)]]))
  matpoints(c(1:4), (test+2^-4), col = 'darkred', type = "p", ylim = lims, log='y', pch=1)
  points(c(1:4), apply(test, 1, mean)+2^-4, col='darkred', type='l', lwd=2.0, lty = lty)
  
  legend("topright", legend = c("untreated", "treated"), col=c('darkblue', "darkred"), bty = "o", lty=c(1, 1))
  
}

dev.off()

##################################################
##################################################
## Section: Enrichment analysis for each stage
##################################################
##################################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

tcs = unique(design$stage)
length(tcs)

source("miRNAseq_functions.R")
list.expressed.miRNAs = read.csv(paste0("../data/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv"), 
                                 header = TRUE, as.is = c(1, 2))

Save.Comparison = TRUE
Check.data.quality = FALSE
Make.new.list.of.expressed.miRNAs = TRUE
Filter.lowly.expressed.using.predefined.miRNA.list = TRUE;
Merge.techinical.replicates.N2 = FALSE

if(Merge.techinical.replicates.N2){
  rep.technical = list(c("57751", "57753"), c("57752", "57754"))
  Merge.techinical.replicates.for.N2(all, id.list=rep.techinical)
}

# make a new list of expressed miRNAs
if(Make.new.list.of.expressed.miRNAs){
  
  kk = which(design$stage == "L1.early")
  design.matrix = data.frame(sample=colnames(all)[-c(1, (kk+1))], design[-kk, ])
  countData = floor(as.matrix(all[, -c(1, (kk+1))]));
  countData[which(is.na(countData))] = 0
  rownames(countData) = all$gene;
  
  #new.list.expressed.miRNA = identify.expressed.miRNAs(countData, design.matrix)  
  new.list.expressed.miRNA = identify.expressed.miRNAs.for.stages(countData, design.matrix)
  #write.csv(new.list.expressed.miRNA, file = paste0(tabDir, "new_list_expressed_miRNAs_lateL1_L4.csv"), row.names=FALSE, quote = FALSE)
}

# filter lowly expressed miRNA with list of predefined miRNAs that were identified using all untreated samples 
if(Filter.lowly.expressed.using.predefined.miRNA.list){
  # prepare old llist
  list.old = find.mature.ones.for.prefixed.expressed.miRNAs(list.expressed.miRNAs)
  #find.mature.ones.for.prefixed.expressed.miRNAs
  # new list
  list.new = find.mature.ones.for.prefixed.expressed.miRNAs(new.list.expressed.miRNA)
  
  old.list = list.old[, c(1:3)]
  new.list = list.new[, c(1:3)]
  old.ggs = unique(old.list$gene)
  new.ggs = unique(new.list$gene)
  cat(setdiff(old.ggs, new.ggs), "\n")
  cat(setdiff(new.ggs, old.ggs), "\n")
  mm = match(old.list$gene, setdiff(old.ggs, new.ggs))
  index2add = which(!is.na(mm))
  expressed.miRNAs = data.frame(rbind(new.list, old.list[index2add, ]), stringsAsFactors = FALSE)
  expressed.miRNAs$mature = expressed.miRNAs$mature>0
  
  save(expressed.miRNAs, file=paste0(RdataDir, 'list_expressed_miRNAs_across_stages', version.analysis, '.Rdata'))
  
}


####################
## start Enrichment Analysis for each neuron class 
####################
require('DESeq2')
read.count = all[, -1];

for(n in 1:length(tcs))
{
  # n = 1
  specifity = tcs[n];
  specDir = paste0(tabDir,specifity, "/")
  if(!dir.exists(specDir)){dir.create(specDir)}
  
  kk = which(design$stage==specifity)
  
  design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
  raw = as.matrix(read.count[,kk])
  raw[which(is.na(raw))] = 0
  raw = floor(raw)
  rownames(raw) = all$gene
  
  ## save plots for enrichment analysis
  pdfname = paste0(specDir, "QCs_assessment_Enrichment_analysis_", specifity, "_", version.analysis, ".pdf") #save all plots during data processing
  pdf(pdfname, width = 12, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  if(Check.data.quality){
    #treat = length(unique(design$treatment[kk]));
    index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
    index.qc = c(1, index.qc, 6)
    
    design.matrix.QC = design.matrix[, index.qc]
    
    if(any(index.qc==3)){
      if(any(design.matrix.QC$genotype=="N2")){
        design.matrix.QC$genotype[which(design.matrix.QC$genotype=="N2")] = "WT"
        if(length(unique(design.matrix.QC$genotype))==1)
          design.matrix.QC = design.matrix.QC[, -which(colnames(design.matrix.QC)=="genotype")]
      }
    }
    
    source("RNAseq_Quality_Controls.R")
    Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix.QC)
  }
  
  ####################
  ## enrichment analysis is done with samples with the same genotype and promoter
  ####################
  # cc = genos[1]; prot = promoters[1]; 
  #cat(specifity, "--", cc, "--", as.character(prot), "\n")
  
  #jj = which(design.matrix$genotype==cc & design.matrix$promoter==prot)
  if(length(unique(design.matrix$treatment))>=2)
  {
    dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ treatment )
    cpm0 = fpm(dds, robust = FALSE)
    
    ## filter lowly expressed miRNA first before estimating scaling factor and dispersion parameters
    if(Filter.lowly.expressed.using.predefined.miRNA.list){
      jj.expressed = NULL
      jj.expressed = match(rownames(dds), expressed.miRNAs$miRNA)
      sels = !is.na(jj.expressed)
      cat("nb of expressed miRNAs --", sum(sels), "\n")
      dds = dds[sels, ]
      index.sel = jj.expressed[sels]
      jj.mature = expressed.miRNAs$mature[index.sel]
      rownames(dds)[jj.mature] = expressed.miRNAs$gene[index.sel[jj.mature]]
      cpm0 = cpm0[sels, ]
      
    }else{sels = c(1:nrow(dds));  dds = dds[sels, ]; } # without filtering
    #cat("size factor is -- ", sizeFactors(dds), "\n")
    
    ## estimate scaling factor and dispersion parameter
    kk.mature = grep("cel-",rownames(dds), invert = TRUE)
    cat('nb of expressed genes -- ', length(kk.mature), "\n")
    dds_sf0 = dds[kk.mature, ]
    dds_sf0 <- estimateSizeFactors(dds_sf0)
    sizeFactors(dds) = sizeFactors(dds_sf0)
    cat("size factor is -- ", sizeFactors(dds), "\n")
    
    dds = estimateDispersions(dds, fitType = "parametric")
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
    
    ## normalization read counts and cpm
    cpm = fpm(dds, robust = TRUE)
    rownames(cpm0) = rownames(cpm);
    
    #ii.test = which(rownames(cpm)=="lsy-6"| rownames(cpm)=="cel-lsy-6-3p");
    #log2(mean(cpm[ii.test, c(3:4)])/mean(cpm[ii.test, c(1:2)]))
    colnames(cpm) = paste0(colnames(cpm), "_",  dds$treatment,  '.normalized.read.counts.DESeq2')
    colnames(cpm0) = paste0(colnames(cpm0), "_",  dds$treatment,  '.cpm')
    dds = nbinomWaldTest(dds, betaPrior = FALSE)
    resultsNames(dds)
    
    res <- results(dds, contrast = c("treatment", "treated", "untreated"));
    summary(res)
    res = as.data.frame(res);
    
    kk.mature = grep("cel-",rownames(res), invert = TRUE)
    plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.8, main=paste0(specifity))
    abline(v=seq(-1, 1, by=0.5), lwd=1.0, col='gray')
    abline(h=c(3, 5, 10), lwd=1.0, lty=2, col='blue')
    text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=1.0, offset = 0.3, pos = 3)
    #ii = match(c('lsy-6', 'mir-791', 'mir-790', 'mir-793'), rownames(res));
    #points(res$log2FoldChange[ii], -log10(res$pvalue)[ii], cex=1.5, col='darkgreen', pch=16)
    
    ## save the comparion in excel form
    if(Save.Comparison){
      #res.sig <- subset(res, padj < FDR.cutoff)
      res.sig = data.frame(cpm0, cpm, res, stringsAsFactors = FALSE);
      #res.sig = res;
      #res.sig = as.data.frame(res.sig);
      res.sig = res.sig[order(-res$log2FoldChange), ]
      #res = res[o1, ]
      #cpm = cpm[o1, ]
      kk.mature = grep("cel-",rownames(res.sig), invert = TRUE)
      kk.star = grep("cel-",rownames(res.sig), invert = FALSE)
      write.csv(res.sig[kk.mature, ], 
                file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, '_Mature_',  version.analysis, '.csv'),     
                row.names = TRUE, quote = FALSE)
      write.csv(res.sig[kk.star, ], 
                file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, '_Star_', version.analysis,'.csv'),     
                row.names = TRUE, quote = FALSE)
    }
  }
  
  dev.off();
  
}

##################################################
##################################################
## Section: Enrichment vs absolute (data normalized with spike-ins)
##################################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_cpm_spikeInsNorm', version.analysis, '.Rdata'))
load(file=paste0(RdataDir, 'list_expressed_miRNAs_across_stages', version.analysis, '.Rdata'))

expressed = expressed.miRNAs[expressed.miRNAs$mature, ]
rownames(expressed) = expressed$miRNA

version.analysis.enrichment = "_miRNAs_timeSeries_spikeIn_R5922_R6016_20180620" 
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis.enrichment, '.Rdata'))
enrich.files = list.files(path = tabDir, pattern = paste0("*Mature_", version.analysis.enrichment, ".csv"), full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

tt = unique(design$stage)
for(n in 1:length(tt))
{
  jj = grep(tt[n], enrich.files)
  test = read.csv(enrich.files[jj], row.names = 1, header = TRUE)
  if(n==1){
    enriched = data.frame(test$pvalue, test$log2FoldChange)
    rownames(enriched) = rownames(test)
  }else{
    mm = match(rownames(enriched), rownames(test))
    enriched = data.frame(enriched, test$pvalue[mm], test$log2FoldChange[mm])
  }
  colnames(enriched)[c(ncol(enriched), (ncol(enriched)-1))] = paste0(tt[n], "_", c("log2FC", "pvalue"))
}
rownames(enriched) = expressed$miRNA[match(rownames(enriched), expressed$gene)]

####################
## Check the global expression in the untreated samples
####################
kk = intersect(grep('normBySpikeIn', colnames(res)), grep('untreated', colnames(res)))
jj = match(rownames(expressed), rownames(res))
xx = res[jj, ]

xx = average.biological.replicates(xx[,kk])
xx = data.frame(log2(xx+10^-4))

pdfname = paste0(resDir, "/heatmap_for_untreatedSamples_normSpikeIns", version.analysis, ".pdf")
pdf(pdfname, width=6, height = 20)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))
library(pheatmap)
library(RColorBrewer)

pheatmap(xx[order(-xx$L1.late_untreated), ], cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

dev.off()

####################
## Check the neuron-specific miRNAs expression
####################
kk = intersect(grep('normBySpikeIn', colnames(res)), grep('_treated', colnames(res)))
xx = res[match(rownames(enriched), rownames(res)), kk]

index.sel = c()
for(n in 1:nrow(enriched))
{
  for(m in 1:length(tt))
  {
    if(m==1){
      if(enriched[n, (2*m-1)] < 0.001 & enriched[n, 2*m] >1) index.sel = c(index.sel, n)
    }else{
      if(enriched[n, (2*m-1)] < 0.01 & enriched[n, 2*m] >1) index.sel = c(index.sel, n)
    }
    
  }
}

index.sel = unique(index.sel)

yy = xx[index.sel, ]
yy = data.frame(average.biological.replicates(yy))

pdfname = paste0(resDir, "/heatmap_for_treatedSamples_normSpikeIns", version.analysis, ".pdf")
pdf(pdfname, width=6, height = 12)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

yy0 = log2(yy)
pheatmap(yy0, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=FALSE, main = "log2(spikeIns.norm)",
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))



max.yy0 = apply(yy0, 1, max)
yy1 = 2^(yy0 - max.yy0)

pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=FALSE, main = "ratio to maximum",
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

dev.off()
#enriched = enriched[, -c(1:2)]
#tt = tt[-1]


#stats.mpu = matrix(NA, nrow = nrow(mpu), ncol=3*length(tt))
#rownames(stats.mpu) = rownames(mpu)
#colnames(stat.mpu) = paste(tt, c("_untreated.mpu", '_treated.mpu', '_ratio.mpu'))


##################################################
##################################################
## Section: session infos
##################################################
##################################################
sessionDir = paste0(resDir, "/log/")
if(!dir.exists(sessionDir)){dir.create(sessionDir)}

sink(paste(sessionDir,"sessionInfo.Chiara.miRNA.DESeq2.txt", sep=""))
sessionInfo()
sink()
save.image(file=paste0(sessionDir,"Chiara_miRNAs_DESeq2.RData"))



