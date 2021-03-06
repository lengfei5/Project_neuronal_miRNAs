##################################################
##################################################
## Project: Chiara's neuron specific miRNA expression
## Script purpose: Identify enriched miRNAs in each neuron class 
## this script is modified on the basis of R script for Chiara et al. 2018 Nature Methods
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Mar  7 10:38:23 2018
##################################################
##################################################
### data verision and analysis version
version.Data = 'miRNAs_neurons_v1';
version.analysis = paste0(version.Data, "_2018_03_07")

miRNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/miRNAseq_functions.R"

### Directories to save results
design.file = "../exp_design/Neuron_project_design_all_v1.xlsx"
dataDir = "../data"
resDir = paste0("../results/miRNAs_neurons_enrichment")
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

Analysis.newData.from.nf_all = TRUE
# Use.N2 = FALSE
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
jj = grep("WT_tax4", design$genotype)
all = process.countTable.v1(all=all, design = design[-jj, ]); # process.countTable.v1 is the function version used before

if(Analysis.newData.from.nf_all){
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
}

xx <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(all, all2))
all = xx


save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

##################################################
##################################################
## Section: Enrichment analysis for each tissue/cell type
##################################################
##################################################
load(file = paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
source("miRNAseq_functions.R")

Filter.lowly.expressed.using.predefined.miRNA.list = FALSE;
Filter.lowly.expressed.using.cpm.threshold = TRUE;
Merge.techinical.replicates.N2 = TRUE

tcs = unique(design$tissue.cell)
tcs = setdiff(tcs, c("whole.body", "whole.body.L3"))
length(tcs)

if(Merge.techinical.replicates.N2){
  rep.technical = list(c("57751", "57753"), c("57752", "57754"), c("66867", "66869"), c("66868", "66870"))
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

####################
## Enrichment Analysis for each neuron class
####################
require('DESeq2')
read.count = all[, -1];

index.N2 = which(design$genotype=="N2")
#design$genotype[index.N2] = "N2"
Save.Comparison = TRUE
Check.data.quality = TRUE

#for(n in 1:length(tcs))
for(n in c(17))
{
  # n = 1
  specifity = tcs[n];
  
  specDir = paste0(tabDir, specifity, "/")
  if(!dir.exists(specDir)){dir.create(specDir)}
  
  kk = which(design$tissue.cell==specifity)
  kk = kk[order(design$promoter[kk])]
  
  if(specifity == "CEPsh.L3"){
    index.N2 = which(design$genotype=="WT" & design$tissue.cell == "whole.body.L3")
  }
  
  # find genome types and promoters
  genos = unique(design$genotype[kk]);
  promoters = unique(design$promoter[kk]); 
  
  if(any(genos=='WT') & specifity != "whole.body") {
    kk = unique(c(kk, index.N2))
  }
  
  design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
  if(specifity == "CEPsh.L3"){
    design.matrix$genotype[which(design.matrix$tissue.cell=="whole.body.L3")] = "N2"
  }
  
  raw = as.matrix(read.count[,kk])
  raw[which(is.na(raw))] = 0
  raw = floor(raw)
  rownames(raw) = all$gene
  
  ## save plots for enrichment analysis
  pdfname = paste0(specDir, "QCs_assessment_Enrichment_analysis_", specifity, ".pdf") #save all plots during data processing
  pdf(pdfname, width = 12, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
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
  if(Check.data.quality){ Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix.QC);}
  
  #dev.off()
  
  ## enrichment analysis is done with samples with the same genotype and promoter
  for(cc in genos)
  {
    for(prot in promoters)
    {
      # cc = genos[1]; prot = promoters[1]; 
      cat(specifity, "--", cc, "--", as.character(prot), "\n")
      
      jj = which(design.matrix$genotype==cc & design.matrix$promoter==prot)
      
      if(length(jj)>=2)
      {
        ## add N2 genotype if WT bakcground for interaction term; 
        if(cc=="WT" & specifity != "whole.body") jj = unique(c(jj, which(design.matrix$genotype=="N2"))) 
        
        countData = raw[, jj]
        if(cc=="WT" & specifity != "whole.body"){
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ genotype + treatment + genotype:treatment)
          dds$genotype <- relevel(dds$genotype, ref="WT");
          dds$treatment = relevel(dds$treatment, ref="untreated");
        }else{
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ treatment )
        }
        
        cpm0 = fpm(dds, robust = FALSE)
        
        ## filter lowly expressed miRNA first before estimating scaling factor and dispersion parameters
        jj.expressed = NULL
        jj.expressed = match(rownames(dds), expressed.miRNAs$miRNA)
        sels = !is.na(jj.expressed)
        cat("nb of expressed miRNAs --", sum(sels), "\n")
        #cat("nb of expressed genes -- ", sum(expressed.miRNAs$mature[sels]))
        dds = dds[sels, ]
        index.sel = jj.expressed[sels]
        jj.mature = expressed.miRNAs$mature[index.sel]
        rownames(dds)[jj.mature] = as.character(expressed.miRNAs$gene[index.sel[jj.mature]])
        cpm0 = cpm0[sels, ]
        
        #cat("size factor is -- ", sizeFactors(dds), "\n")
        
        ## estimate scaling factor and dispersion parameter
        kk.mature = grep("cel-",rownames(dds), invert = TRUE)
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
        
        if(cc=="WT" & specifity != "whole.body")
        {
          ## the ones showing difference treated vs untreated in WT background
          res = results(dds, contrast=c("treatment","treated","untreated"))
          summary(res)
          res0 = res;
          colnames(res0) = paste0(colnames(res0), ".without.N2")
          kk.mature = grep("cel-",rownames(res), invert = TRUE)
          plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.8, 
               main=paste0(cc, "--", prot, " (NOT using N2)"))
          abline(v=0, lwd=2.0, col='black')
          abline(h=c(5, 10), lwd=1.0, lty=2, col='blue')
          text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=0.7, offset = 0.3, pos = 3)
          
          ## the one showing N2-specific treatment effect is smaller than WT-specific treatment effect 
          ## (refer to the DEseq2 mannul or tutorial https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
          res2 = results(dds, name = "genotypeN2.treatmenttreated", lfcThreshold = 0, altHypothesis = "less")
          
          ## merge both res1 and res2
          res1 = as.data.frame(res); 
          res2 = as.data.frame(res2);
          
          ## change the sign (log2(FC.WT.treated) - log2(FC.WT.untreated)) - (log2(FC.N2.treated) - log2(FC.N2.untreated))
          res2$log2FoldChange = -res2$log2FoldChange; 
          ii.enrich = which(res1$log2FoldChange>0) 
          res1[ii.enrich, ] = res2[ii.enrich, ] ## replaced by res2 if enriched; keep if depleted
          
          ## merge the table with comparison with N2 and the one without
          res = data.frame(res1, res0, stringsAsFactors = FALSE)
          
        }else{
          res <- results(dds, contrast = c("treatment", "treated", "untreated"));
          summary(res)
          res = as.data.frame(res);
        }
        
        kk.mature = grep("cel-",rownames(res), invert = TRUE)
        plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', 
             cex=0.8, main=paste0(cc, "--", prot))
        abline(v=seq(-1, 1, by=0.5), lwd=1.0, col='gray')
        abline(h=c(3, 5, 10), lwd=1.0, lty=2, col='blue')
        text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=0.7, offset = 0.3, pos = 3)
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
                    file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Mature_',  version.analysis, '.csv'),     
                    row.names = TRUE, quote = FALSE)
          write.csv(res.sig[kk.star, ], 
                    file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Star_', version.analysis,'.csv'),     
                    row.names = TRUE, quote = FALSE)
        }
      }
    }
  
  }
  
  dev.off();
  
}
  
###############################
# log session info
###############################
sessionDir = paste0(resDir, "/log/")
if(!dir.exists(sessionDir)){dir.create(sessionDir)}

sink(paste(sessionDir,"sessionInfo.Chiara.miRNA.DESeq2", version.analysis, ".txt", sep=""))
sessionInfo()
sink()
save.image(file=paste0(sessionDir,"Chiara_miRNAs_DESeq2", version.analysis, ".RData"))

########################################################
########################################################
# Section: other additional plots or test
########################################################
########################################################
Plot.for.presentation = FALSE
if(Plot.for.presentation){
  yy = xx[grep('mir-789-2', rownames(xx)), ]
  yy = yy[,grep('normalized.read.counts.DESeq2', colnames(yy))]
  yy = unlist(yy)
  for(n in 1:5)
  {
    cols = 'darkblue'
    lwd = 2.0
    if(n==1) plot(c(1:2), yy[c((2*n-1), 2*n)], ylim=range(yy), type= 'b', log='y', col=cols, lwd=3.0)
    else {
      if(n >2) {cols = 'red'; lwd=1.0}
      points(c(1:2), yy[c((2*n-1), 2*n)], type='b', lwd=lwd, col= cols)
    }
  }
}

###############################
# compare the wt and henn1.mutant and the relationship    
###############################
Compare.rab3_wt_vs_henn1.mutant = FALSE
if(Compare.rab3_wt_vs_henn1.mutant){
  test.Dir = paste0(resDir, "/tables/ASE.neurons")
  xlist = list.files(path = test.Dir, pattern = "*.csv", full.names = TRUE)
  xlist = xlist[grep("_Mature", xlist)]
  
  xx = read.csv(xlist[grep("henn1.mutant", xlist)], header = TRUE, row.names = 1)
  jj1 = intersect(grep("DESeq2", colnames(xx)), grep("_treated", colnames(xx)))
  jj2 = intersect(grep("DESeq2", colnames(xx)), grep("_untreated", colnames(xx)))
  plot(log2(apply(as.matrix(xx[, jj1]), 1, mean)/apply(as.matrix(xx[, jj2]), 1, mean)), xx$log2FoldChange, 
       xlab = "log2FC from my own calculation", ylab = "log2FC from DESeq2", main = basename(test.Dir))
  abline(0, 1, lwd=2.0, col='red')
  
  yy = read.csv(xlist[grep("WT", xlist)], header = TRUE, row.names = 1)
  
  mm = match(rownames(xx), rownames(yy))
  yy = yy[mm, ]
  
  source("miRNAseq_functions.R")
  
  kk = intersect(grep("DESeq2", colnames(yy)), grep("N2", colnames(yy)))
  fc.N2 = log2(apply(yy[, kk[grep("_treated_", colnames(yy)[kk])]], 1, mean)/
                 apply(yy[, kk[grep("_untreated_", colnames(yy)[kk])]], 1, mean))
  
  kk = intersect(grep("DESeq2", colnames(yy)), grep("WT", colnames(yy)))
  fc.wt = log2(apply(yy[, kk[grep("_treated_", colnames(yy)[kk])]], 1, mean)/
                 apply(yy[, kk[grep("_untreated_", colnames(yy)[kk])]], 1, mean))
  
  par(mfrow=c(2, 2))
  plot(xx$log2FoldChange, fc.wt, xlab='henn1.mutant', ylab = "WT")
  abline(0, 1, lwd=2.0, col='red')
  #text( xx$log2FoldChange, fc.wt, rownames(xx), offset = 0.7, cex = 0.5)
  
  #kk = c(1:40)
  #plot(xx$log2FoldChange, yy$log2FoldChange.without.N2)
  #abline(0, 1, lwd=2.0, col='red')
  
  plot(xx$log2FoldChange, yy$log2FoldChange, ylim = c(-4, 10))
  abline(0, 1, lwd=2.0, col='red')
  
  plot(fc.N2, fc.wt, xlab='N2', ylab = "WT")
  abline(0, 1, lwd=2.0, col='red')
  
  #jj = which(xx$log2FoldChange- fc.wt<0)
  #plot()
  
  plot(xx$log2FoldChange, (fc.wt-fc.N2), xlab='henn1.mutant', ylab = "WT-N2")
  abline(0, 1, lwd=2.0, col='red')
  #text( xx$log2FoldChange, fc.wt, rownames(xx), offset = 0.7, cex = 0.5)
  
}
