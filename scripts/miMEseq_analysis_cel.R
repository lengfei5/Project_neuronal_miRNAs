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

#### Functions
process.countTable = function(all, design)
{
  index = c()
  for(n in 1:nrow(design))
  {
    #n = 1;
    jj = intersect(grep(design$SampleID[n], colnames(all)), grep("Total.count", colnames(all)));
    if(length(jj)==1) {
      index = c(index,jj)
    }else{print(paste0("ERROR for sample--", design$SampleID[n]))}
  }
  
  newall = data.frame(as.character(all[,1]),  as.matrix(all[, index]), stringsAsFactors = FALSE)
  colnames(newall)[1] = "gene";
  colnames(newall)[-1] = paste0(design$genotype, "_", design$tissue.cell, "_", design$SampleID)
  
  return(newall)
}

### data verision and analysis version   
version.Data = 'miRNAs_cel';
version.analysis = paste0("_", version.Data)

### Directories to save results 
dataDir = "../data/"
resDir = "../results/cel_cpm_normalizedCount/"
tabDir =  paste0(resDir, "tables/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

#### Import Sample information and table of read counts
library("openxlsx")
design = read.xlsx(paste0(dataDir, "Sample_Information_cel.xlsx"), sheet = 1, colNames = TRUE)
design$genotype[grep("WT", design$genotype)] = "WT"
design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole_animal_no_promoter"

all = read.delim(paste0(dataDir, "cel_countTable.txt"), sep = "\t", header = TRUE)
# processing count table
all = process.countTable(all=all, design = design);

#### Enrichment analysis for each tissue/cell type
tcs = unique(design$tissue.cell)
tcs = setdiff(tcs, c("whole_animal", "all_tissues"))
length(tcs)

# filter lowly expressed miRNA with list of predefined miRNAs that were identified using all untreated samples 
Filter.lowly.expressed.using.predefined.miRNA.list = TRUE;
if(Filter.lowly.expressed.using.predefined.miRNA.list){
  list.expressed = read.csv(paste0(dataDir, "list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv"), header = TRUE, as.is = c(1, 2))
  colnames(list.expressed)[c(1,2)] = c("miRNA", "gene")
  list.expressed = list.expressed[order(list.expressed$miRNA), ]
  ggs.uniq = unique(list.expressed$gene)
  
  mature = rep(NA, nrow(list.expressed))
  kk = grep('.cpm', colnames(list.expressed))
  for(n in 1:length(ggs.uniq))
  {
    jj = which(list.expressed$gene==ggs.uniq[n])
    if(length(jj)>1){
      index.max = apply(list.expressed[jj, kk], 2, which.max)
      nb.first.max = length(which(index.max==1))
      nb.sec.max = length(which(index.max==2))
      #cat(n, ": ", as.character(ggs.uniq[n]), "--",  nb.first.max, "--", nb.sec.max, "\n")
      if(nb.first.max>nb.sec.max){
        mature[jj[1]] = TRUE; mature[jj[2]] = FALSE;  
      }else{
        mature[jj[1]] = FALSE; mature[jj[2]] = TRUE; 
      }
    }else{
      #cat(n,": ", as.character(ggs.uniq[n]),  "-- no selection \t")
      mature[jj] = TRUE;
    }
  }
  expressed.miRNAs = data.frame((list.expressed[, c(1, 2)]), mature=(mature), list.expressed[, -c(1:2)], stringsAsFactors = FALSE)
}

read.count = all[, -1];
require('DESeq2')
#index.N2 = which(design$genotype=="N2")
index.N2 = grep("N2", design$sampleInfo)
design$genotype[index.N2] = "N2"
Save.Comparison = TRUE


for(n in 1:length(tcs))
{
  # n = 1
  specifity = tcs[n];
  kk = which(design$tissue.cell==specifity)
  
  # find genome types and promoters
  genos = unique(design$genotype[kk]);
  promoters = unique(design$promoter[kk]); 
  
  if(any(genos=='WT') & specifity != "whole_animal")  
  {
    kk = unique(c(kk, index.N2))
  }
  
  design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
  raw = as.matrix(read.count[,kk])
  raw[which(is.na(raw))] = 0
  raw = floor(raw)
  rownames(raw) = all$gene
  
  ## save plots for enrichment analysis
  pdfname = paste0(resDir, "Enrichment_analysis_", specifity, ".pdf") #save all plots during data processing
  pdf(pdfname, width = 12, height = 8)
  
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
        if(cc=="WT" & specifity != "whole_animal") jj = unique(c(jj, which(design.matrix$genotype=="N2"))) 
        
        countData = raw[, jj]
        if(cc=="WT" & specifity != "whole_animal"){
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ genotype + treatment + genotype:treatment)
          dds$genotype <- relevel(dds$genotype, ref="WT");
          dds$treatment = relevel(dds$treatment, ref="untreated");
        }else{
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ treatment )
        }
        
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
        dds_sf0 = dds[kk.mature, ]
        dds_sf0 <- estimateSizeFactors(dds_sf0)
        sizeFactors(dds) = sizeFactors(dds_sf0)
        cat("size factor is -- ", sizeFactors(dds), "\n")
        
        dds = estimateDispersions(dds, fitType = "parametric")
        plotDispEsts(dds, ylim=c(0.001, 10))
        
        ## normalization read counts and cpm
        cpm = fpm(dds, robust = TRUE)
        rownames(cpm0) = rownames(cpm);
        
        #ii.test = which(rownames(cpm)=="lsy-6"| rownames(cpm)=="cel-lsy-6-3p");
        #log2(mean(cpm[ii.test, c(3:4)])/mean(cpm[ii.test, c(1:2)]))
        colnames(cpm) = paste0(colnames(cpm), "_",  dds$treatment,  '.normalized.read.counts.DESeq2')
        colnames(cpm0) = paste0(colnames(cpm0), "_",  dds$treatment,  '.cpm')
        dds = nbinomWaldTest(dds, betaPrior = FALSE)
        resultsNames(dds)
        
        if(cc=="WT" & specifity != "whole_animal")
        {
          ## the ones showing difference treated vs untreated in WT background
          res = results(dds, contrast=c("treatment","treated","untreated"))
          summary(res)
          res0 = res;
          colnames(res0) = paste0(colnames(res0), ".without.N2")
          kk.mature = grep("cel-",rownames(res), invert = TRUE)
          plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.7, main=paste0(cc, "--", prot, " (NOT using N2)"))
          abline(v=0, lwd=2.0, col='black')
          abline(h=c(5, 10), lwd=1.0, lty=2, col='blue')
          text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=0.7, offset = 0.3, pos = 3)
          
          ## the one showing N2-specific treatment effect is smaller than WT-specific treatment effect 
          ## (refer to the DEseq2 mannul or tutorial https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
          res2 = results(dds, name = "genotypeN2.treatmenttreated", lfcThreshold = 0, altHypothesis = "less")
          
          ## merge both res1 and res2
          res1 = as.data.frame(res); 
          res2 = as.data.frame(res2);
          res2$log2FoldChange = -res2$log2FoldChange; ## change the sign (log2(FC.WT.treated) - log2(FC.WT.untreated)) - (log2(FC.N2.treated) - log2(FC.N2.untreated))
          ii.enrich = which(res1$log2FoldChange>0) 
          res1[ii.enrich, ] = res2[ii.enrich, ] ## replaced by res2 if enriched; keep if depleted
          
          ## merge the table with comparison with N2 and the one without
          #res = data.frame(res1, res0, stringsAsFactors = FALSE)
        }else{
          res <- results(dds, contrast = c("treatment", "treated", "untreated"));
          summary(res)
          res = as.data.frame(res);
        }
        
        kk.mature = grep("cel-",rownames(res), invert = TRUE)
        plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.7, main=paste0(cc, "--", prot))
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
                    file=paste0(tabDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Mature',  version.analysis, '.csv'),     
                    col.names = TRUE, row.names = TRUE, quote = FALSE)
          write.csv(res.sig[kk.star, ], 
                    file=paste0(tabDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Star', version.analysis,'.csv'),     
                    col.names = TRUE, row.names = TRUE, quote = FALSE)
        }
      }
    }
  
  }
  
  dev.off();
}

####################
## log session info
####################
sessionDir = paste0(resDir, "log/")
if(!dir.exists(sessionDir)){dir.create(sessionDir)}

sink(paste(sessionDir,"sessionInfo.Chiara.miRNA.DESeq2.txt", sep=""))
sessionInfo()
sink()
save.image(file=paste0(sessionDir,"Chiara_miRNAs_DESeq2.RData"))
