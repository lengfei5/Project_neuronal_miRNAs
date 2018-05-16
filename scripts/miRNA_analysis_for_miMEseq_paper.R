########################################################
#### verision of data and analysis    
########################################################
version.Data = 'miRNAs_cel';
version.analysis = paste0("2017-11-15_", version.Data)

dataDir = "../final/data/"
resDir = "../final/restuls/cel/"
if(!dir.exists(resDir)){dir.create(resDir)}

########################################################
########################################################
#### Import data and sample information    
########################################################
########################################################
library("openxlsx")

xx = read.xlsx("../final/HEN1_project_design_all_4paper.xlsx", sheet = 1, colNames = TRUE)

design = data.frame(xx$ID, as.character(xx$genetic.background), xx$tissue.cell.type, xx$promoters, xx$treatment, xx$Sample, stringsAsFactors = FALSE)
colnames(design) = c('SampleID', 'genotype', 'tissue.cell', 'promoter', 'treatment', 'sampleInfo')
kk = grep("No Treatment", design$treatment)
design$treatment[kk] = "untreated"
kk = grep("Treatment", design$treatment)
design$treatment[kk] = "treated"
jj = grep("henn-1", design$genotype)
design$genotype[jj] = "henn1_mutant"

write.xlsx(design, file=paste0(dataDir, "Sample_Information_cel.xlsx"))


##### Import tables
#data.path = paste0('DATA')
xlist<-list.files(path=paste0('../DATA'), pattern = "*.CountAndUMIfrRaw.txt", full.names = TRUE)
ylist = list.files(path=paste0('../DATA'), pattern = "*.typeHierarchy.txt", full.names = TRUE)

counts = list()
for(n in 1:length(xlist)){
  ff = read.table(xlist[n], sep='\t', header = TRUE);
  counts[[n]] = data.frame(ff, stringsAsFactors = FALSE);
};
stats = NULL
for(n in 1:length(ylist)){
  #n = 1;
  xx = read.table(ylist[n], sep='\t', header = TRUE);
  xx = xx[which(!is.na(xx[,1])==TRUE), ]
  xx = data.frame(xx, stringsAsFactors = FALSE);
  if(n == 1){
    stats = xx;
  }else{
    mm = match(stats$type, xx$type)
    stats = data.frame(stats, xx[mm, -1])
  }
  #stats[[n]] = data.frame(xx, stringsAsFactors = FALSE);
};

### make the data table
find.mirName = function(x){test = unlist(strsplit(as.character(x), '-'));
  return(paste0(test[-c(1,length(test))], collapse = '-'))
} # this function is not used here

all = NULL
ggs = NULL
norms = NULL
Select.Abundant.Arms = FALSE
for(n in 1:nrow(design))
{
  #cat(design[n, ], "\n")
  if(n == 1){
    for(m in 1:length(counts))
    {
      data = counts[[m]]
      kk = grep(design$SampleID[n], colnames(data))
      if(length(kk)>0){
        cat(m, "\n")
        data = data[, c(1:5, kk)]
        
        ## previous version analysis was to select abundant arm with the following code, but new version will keep all arms
        if(Select.Abundant.Arms){
          jj = which(colnames(data)==paste0('miR.count.', design$SampleID[n]));
          data = data[which(data[,jj] & !is.na(data[,jj])),] 
          ggs = sapply(data$gname, find.mirName);
        }else{ggs = data$gname}
        all = data.frame(ggs, data[, c(grep('count.total', colnames(data)), grep('UMIfr.total', colnames(data)))], stringsAsFactors = FALSE)
        #all = data.frame(ggs, data[, grep('count.total', colnames(data))], data[, grep('UMIfr.total', colnames(data))], stringsAsFactors = FALSE)
      }
    }
    norms = data.frame(stats$type, stats[, grep(design$SampleID[n], colnames(stats))], stringsAsFactors = FALSE)
    colnames(norms)[1] = "type"
    colnames(norms)[(n+1)] = design$SampleID[n]
    colnames(all)[1] = 'gene'
  }else{
    for(m in 1:length(counts))
    {
      data = counts[[m]]
      kk = grep(design$SampleID[n], colnames(data))
      if(length(kk)>0){
        data = data[, c(1:5, kk)]
        Select.Abundant.Arms = FALSE
        if(Select.Abundant.Arms){
          jj = grep('miR.count', colnames(data))
          data = data[which(data[,jj] & !is.na(data[,jj])),]
          ggs = sapply(data$gname, find.mirName);
        }else{
          ggs = data$gname;
        }
        gg.union = unique(union(all$gene, ggs))
        all = data.frame(gg.union, 
                         all[match(gg.union, all$gene), -1],
                         data[match(gg.union, ggs), c(grep('count.total', colnames(data)), grep('UMIfr.total', colnames(data)))], 
                         stringsAsFactors = FALSE)
      }
    }
    
    norms = data.frame(norms, stats[, grep(design$SampleID[n], colnames(stats))], stringsAsFactors = FALSE)
    #colnames(norms)[1] = "type"
    colnames(all)[1] = 'gene'
    colnames(norms)[(n+1)] = design$SampleID[n]
    print(all$gene[which(is.na(all$gene))])
  }
}

colnames(norms)[-1] = paste0('readcounts_', design$SampleID)
jj = which(c(2:ncol(all))%%2==0)+2
UMI = all[, c(1, jj)]
all = all[, c(1, (jj-1))]
colnames(all)[-1] = paste0(design$genotype, "_", design$tissue.cell, "_", design$promoter, "_", design$SampleID)
colnames(UMI)[(-1)] = paste0(design$genotype, "_", design$tissue.cell, "_", design$promoter, "_", design$SampleID, "_UMI")

jj = grep('lsy-6', all$gene)
print(all[jj, ])

save(design, all, norms, UMI, file=paste0('../Rdata/', version.Data, '_Design_Raw_readCounts_UMI_norms.Rdata'))

########################################################
########################################################
#### Check quality for each tissue
########################################################
########################################################
load(file=paste0('../Rdata/', version.Data, '_Design_Raw_readCounts_UMI_norms.Rdata'))

### Check relationship between read counts and UMI
Check.UMI.Read.Counts = TRUE
if(Check.UMI.Read.Counts){
  pdfname = paste0("PLOTs/", version.analysis, "_Read_Counts_UMI_relationship.pdf")
  pdf(pdfname, width = 10, height =8 )
  
  for(n in 2:ncol(all))
  {
    plot(UMI[,n], all[,n], log='xy', xlab='nb of UMI', ylab='raw read counts', main=colnames(all)[n])
    abline(0, 1, lwd=2.0, col='blue')
  }
  dev.off()
}

########################################################
########################################################
#### control the quality and enrichment analysis for each tissue/cell type
########################################################
########################################################
## import data
load(file=paste0('../Rdata/', version.Data, '_Design_Raw_readCounts_UMI_norms.Rdata'))
# design = design[which(design$SampleID!="36229" & design$SampleID !="36211"), ]

## analysis version and folder for the results
#version.analysis = paste0("2017-09-15_", version.Data)
#resDir = "../Results/2017-09-15/"
#tabDir =  paste0(resDir, "Tables/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

## check the fraction of miRNAs which change the arms after treatment
Check.3p.5p.arm.switch = FALSE
if(Check.3p.5p.arm.switch) Check.3p.5p.arm.switch(all);

length(unique(design$tissue.cell))
tcs = unique(design$tissue.cell)
tcs = setdiff(tcs, c("Glial_cells"))

require('DESeq2')
index.N2 = grep('N2', design$sampleInfo);

Save.Comparison = TRUE

Filter.lowly.expressed.using.predefined.miRNA.list = TRUE;
if(Filter.lowly.expressed.using.predefined.miRNA.list){
  #list.expressed = read.csv("../DATA/list_expressed_miRNAs_using_cpm_Untreated_samples_Henn1_mutant_WT_all_by_Luisa_Chiara.csv", header = TRUE)
  list.expressed = read.csv("../DATA/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv", header = TRUE, as.is = c(1, 2))
  colnames(list.expressed)[c(1,2)] = c("miRNA", "gene")
  list.expressed = list.expressed[order(list.expressed$miRNA), ]
  ggs.uniq = unique(list.expressed$gene)
  
  mature = rep(NA, nrow(list.expressed))
  kk = grep('.cpm', colnames(list.expressed))
  for(n in 1:length(ggs.uniq))
  {
    #n = 1
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
  #list.expressed = expre$miRNA;
  #list.expressed = list.expressed[which(list.expressed != "" & !is.na(list.expressed))]
}

Save.all.nontreated.cpm = FALSE;
if(Save.all.nontreated.cpm){all.untreated = NULL; names.untreated=NULL;}
#filter.lowly.expressed..with.bimodal.distribution = FALSE
#Compare.WT.henn1.mutant = FALSE

read.count = all[, -1];
source("RNAseq_Quality_Controls.R")

for(n in 1:6)
#for(n in 1:length(tcs))
{
  # n = 2
  specifity = tcs[n];
  kk = which(design$tissue.cell==specifity)
  norms_scaled = norms[, c(1, (1+kk))]; 
  
  genos = unique(design$genotype[kk]);
  promoters = unique(design$promoter[kk]); 
  
  if(any(genos=='WT') & specifity != "whole_animal")  kk = unique(c(kk, index.N2))
  
  design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
  raw = as.matrix(read.count[,kk])
  raw[which(is.na(raw))] = 0
  ### start enrichment analysis 
  raw = floor(raw)
  rownames(raw) = all$gene
  
  #treat = length(unique(design$treatment[kk]));
  index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  index.qc = c(1, index.qc, 6)
  
  source("RNAseq_Quality_Controls.R")
  pdfname = paste0(resDir, "Data_qulity_assessment_", specifity, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix[, index.qc])
  dev.off()
  
  pdfname = paste0(resDir, "Enrichment_analysis_", specifity, ".pdf") #save all plots during data processing
  pdf(pdfname, width = 12, height = 8)
  
  ## loop for genos and promoters
  ## DE analysis is done with samples with the same genotype and promoter
  jj.N2 = NULL;
  if(specifity != "whole_animal")
  {
    jj.N2 =  grep("N2", design.matrix$sampleInfo)
    if(length(jj.N2)>0) design.matrix$genotype[jj.N2] = "N2"
  }
  
  for(cc in genos)
  {
    for(prot in promoters)
    {
      # cc = genos[1]; prot = promoters[1]; 
      cat(specifity, "--", cc, "--", as.character(prot), "\n")
      jj = which(design.matrix$genotype==cc & design.matrix$promoter==prot)
      if(length(jj)>=2)
      {
        ## discard bad samples for phyrynx and rps-5 promoter in henn1-mutant
        if(specifity=="Pharynx" & cc=="henn1_mutant" & prot=="myo-2")  ## remove bad samples from 1st sequencing for Pharynx
        {
          jj = setdiff(jj, match(c("34094", "34095", "34096", "34097"), design.matrix$SampleID));
        }
        if(specifity=="whole_animal" & cc=="henn1_mutant" & prot=="rps-5") ## remove bad samples for rps-5 promter in henn1_mutant background
        {
          jj = setdiff(jj, match(c("36229", "36211"), design.matrix$SampleID)); 
          # henn1_mutant             whole_animal       rps-5
          #36229 and 36211
        }
        
        ## add N2 genotype in the WT background; 
        if(cc=="WT" & specifity != "whole_animal") jj = unique(c(jj, jj.N2)) 
        
        countData = raw[, jj]
        if(cc=="WT" & specifity != "whole_animal"){
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ genotype + treatment + genotype:treatment)
          dds$genotype <- relevel(dds$genotype, ref="WT");
          dds$treatment = relevel(dds$treatment, ref="untreated");
        }else{
          dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix[jj, ]), design = ~ treatment )
        }
        
        ####
        #### CONSIDERATION here : normlization the data before or after filtering lowly expressed genes
        #### this is extremely tricky, because DESeq2 normalization is affected to the filtering and secondly the hypothesis here is not applicable.
        #### After testing, just normalization using only expressed and mature arms of miRNAs may work the best
        ####
        #dds0 = dds
        #dds1 <- estimateSizeFactors(dds0)
        #cat("size factor is -- ", sizeFactors(dds1), "\n")
        
        ## filter lowly expressed miRNA first before estimating scaling factor and dispersion parameters
        if(Filter.lowly.expressed.using.predefined.miRNA.list){
          ## here use the curated list of expressed miRNAs by Luisa and Chiara
          jj.expressed = NULL
          jj.expressed = match(rownames(dds), expressed.miRNAs$miRNA)
          sels = !is.na(jj.expressed)
          cat("nb of expressed miRNAs --", sum(sels), "\n")
          dds = dds[sels, ]
          index.sel = jj.expressed[sels]
          jj.mature = expressed.miRNAs$mature[index.sel]
          rownames(dds)[jj.mature] = expressed.miRNAs$gene[index.sel[jj.mature]]
          
        }else{sels = c(1:nrow(dds));  dds = dds[sels, ]; } # without filtering
        #cat("size factor is -- ", sizeFactors(dds), "\n")
        
        ## estimate scaling factor and dispersion parameter
        #dds2 <- estimateSizeFactors(dds)
        #cat("size factor is -- ", sizeFactors(dds2), "\n")
        kk.mature = grep("cel-",rownames(dds), invert = TRUE)
        dds_sf0 = dds[kk.mature, ]
        dds_sf0 <- estimateSizeFactors(dds_sf0)
        sizeFactors(dds) = sizeFactors(dds_sf0)
        cat("size factor is -- ", sizeFactors(dds), "\n")
        
        dds = estimateDispersions(dds, fitType = "parametric")
        plotDispEsts(dds, ylim=c(0.001, 10))
        
        cpm = fpm(dds, robust = TRUE)
        #ii.test = which(rownames(cpm)=="lsy-6"| rownames(cpm)=="cel-lsy-6-3p");
        #log2(mean(cpm[ii.test, c(3:4)])/mean(cpm[ii.test, c(1:2)]))
        colnames(cpm) = paste0(colnames(cpm), "_",  dds$treatment,  '.cpm')
        
        if(Save.all.nontreated.cpm){new.untreated = setdiff(colnames(cpm)[which(dds$treatment=="untreated")], names.untreated);
          if(length(new.untreated)>0){ii.new = match(new.untreated, colnames(cpm));
            all.untreated = cbind(all.untreated, cpm[, ii.new]);
            names.untreated = c(names.untreated, new.untreated); } }
        
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
          ## log2(FC.N2.treated) - log2(FC.N2.untreated) - (log2(FC.WT.treated) - log2(FC.WT.untreated))
          res2 = results(dds, name = "genotypeN2.treatmenttreated", lfcThreshold = 0, altHypothesis = "less")
          
          ## merge both res1 and res2
          res1 = as.data.frame(res); 
          res2 = as.data.frame(res2);
          res2$log2FoldChange = -res2$log2FoldChange; ## change the sign (log2(FC.WT.treated) - log2(FC.WT.untreated)) - (log2(FC.N2.treated) - log2(FC.N2.untreated))
          ii.enrich = which(res1$log2FoldChange>0) 
          res1[ii.enrich, ] = res2[ii.enrich, ] ## replaced by res2 if enriched; keep if depleted
          
          ## merge the table with comparison with N2 and the one without
          res = data.frame(res1, res0, stringsAsFactors = FALSE)
          # genes which have a different genotype effect in time 1 than in time 0
          #res.KOvsWT.int <- results(dds2, name = "genotypeko.time400min", lfcThreshold = LFC.cutoff, altHypothesis = "less")
          #colnames(res.KOvsWT.int) = paste0(colnames(res.KOvsWT.int), ".400vs90min.WT.increase.more")
        }else{
          res <- results(dds, contrast = c("treatment", "treated", "untreated"));
          #res <- results(dds, pAdjustMethod='BH', alpha = FDR.cutoff, lfcThreshold = LFC.cutoff, contrast = c("conds", cc))
          #table(res$padj < FDR.cutoff)
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
          res.sig = data.frame(cpm, res, stringsAsFactors = FALSE);
          #res.sig = res;
          #res.sig = as.data.frame(res.sig);
          res.sig = res.sig[order(-res$log2FoldChange), ]
          #res = res[o1, ]
          #cpm = cpm[o1, ]
          kk.mature = grep("cel-",rownames(res.sig), invert = TRUE)
          kk.star = grep("cel-",rownames(res.sig), invert = FALSE)
          write.csv(res.sig[kk.mature, ], 
                    file=paste0(tabDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Mature.csv'),     
                    col.names = TRUE, row.names = TRUE, quote = FALSE)
          write.csv(res.sig[kk.star, ], 
                    file=paste0(tabDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Star.csv'),     
                    col.names = TRUE, row.names = TRUE, quote = FALSE)
        }
      }
    }
  
  }
  
  dev.off();
  
}

if(Save.all.nontreated.cpm){
  colnames(all.untreated) = names.untreated
  
  rr = matrix(NA, ncol=ncol(all.untreated), nrow = nrow(all.untreated));
  for(n in 1:ncol(rr)) rr[,n] = order(all.untreated[,n], decreasing = TRUE);
  #rr = apply(all.untreated, 2, order)
  #rr = cbind(rr, apply(rr, 1, mean))
  xx = apply(all.untreated, 1, mean)
  
  xx = data.frame(all.untreated, xx)
  
  colnames(xx)[ncol(xx)] = 'mean'
  xx = xx[order(-xx$mean), ]
  
  ggs = sapply(rownames(xx), find.mirName);
  
  xx = data.frame(ggs, xx, stringsAsFactors = FALSE)
  write.csv(xx, file=paste0(tabDir, 'miRNA_all_untreated_3p_5p.csv'),     
            col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  ####
  ### 10 cpm as threshold
  jj = which(xx$mean>10)
  length(unique(xx$ggs[jj]))
  ggs.uniq = unique(xx$ggs[jj]);
  
  setdiff(list.expressed, ggs.uniq)
  setdiff(ggs.uniq, list.expressed)
  mm = match(list.expressed, unique(xx$ggs))
  
  mm = match(xx$ggs, ggs.uniq)
  yy = xx[which(!is.na(mm)), ]
  
  write.csv(yy, file="../DATA/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

########
## log session info
########
sessionDir = paste0(resDir, "log/")
if(!dir.exists(sessionDir)){dir.create(sessionDir)}

sink(paste(sessionDir,"sessionInfo.Chiara.miRNA.DESeq2.txt", sep=""))
sessionInfo()
sink()
save.image(file=paste0(sessionDir,"Chiara_miRNAs_DESeq2.RData"))
