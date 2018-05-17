##########################################################################
##########################################################################
## Project: Chiara's neuron class-specific miRNA expression  
## Script purpose: to test piRNA normalization for the expression matrix in the decomvolution analysis
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed May 16 13:58:52 2018
##########################################################################
##########################################################################
RdataDir = paste0("../results/miRNAs_neurons_v1_2018_03_07/Rdata/")
statDir = "../data/normalized_piRNAs"
version.table = "miRNAs_neurons_v1_2018_03_07"
resDir = "../results/tables_for_decomvolution"
if(!dir.exists(resDir)) dir.create(resDir)
Save.Processed.Tables = TRUE
######################################
######################################
## Section: load tables and mapping statistics for piRNA and siRNAs
######################################
######################################
load(file = paste0(RdataDir, 'Design_Raw_readCounts_', version.table, '.Rdata'))
source("miRNAseq_functions.R")

Filter.lowly.expressed.using.predefined.miRNA.list = TRUE;
Merge.techinical.replicates.N2 = TRUE

tcs = unique(design$tissue.cell)
tcs = setdiff(tcs, c("whole.body"))
length(tcs)

if(Merge.techinical.replicates.N2){
  rep.technical = list(c("57751", "57753"), c("57752", "57754"))
  for(n in 1:length(rep.technical))
  {
    index = c()
    for(id in rep.technical[[n]])
    {
      #print(id)
      index = c(index, which(design$SampleID==id))
    }
    
    design$SampleID[index[1]] = paste0(design$SampleID[index], collapse = ".")
    ss = apply(all[, (index+1)], 1, function(x) sum(x, na.rm = TRUE))
    all[, (index[1]+1)] = ss;
    colnames(all)[(index[1]+1)] = paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  design$SampleID[index[1]])
    design = design[-index[-1], ]
    all = all[, -(index[-1]+1)]
  }
} 

# filter lowly expressed miRNA with list of predefined miRNAs that were identified using all untreated samples 
if(Filter.lowly.expressed.using.predefined.miRNA.list){
  list.expressed = read.csv(paste0("../data/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv"), header = TRUE, as.is = c(1, 2))
  # prepare old llist
  list.expressed = find.mature.ones.for.expressed.miRNAs(list.expressed)
  expressed.miRNAs = data.frame(list.expressed[, c(1:3)], stringsAsFactors = FALSE)
}

####################
## calculate cpm and select only expressed miRNAs and start to use the gene names instead of arms 
####################
raw = as.matrix(all[, -1])
raw[which(is.na(raw))] = 0
raw = floor(raw)
rownames(raw) = all$gene
library.sizes = apply(raw, 2, sum)

#cpm = my.cpm.normalization(raw)

expressed.miRNAs = expressed.miRNAs[expressed.miRNAs$mature,]
mm = match(expressed.miRNAs$miRNA, all$gene)
countData = raw[mm, ]
rownames(countData) = expressed.miRNAs$gene

## filter the samples unrelevant (non-neuron samples, henn1-mutant background)
kk = which(design$tissue.cell=="Glial.cells" | design$tissue.cell == "CEPsh" | design$genotype == "henn1.mutant")
if(length(kk)>0){
  design.matrix = design[-kk, ]
  countData = countData[, -kk]
  library.sizes = library.sizes[-kk]
}

######################################
######################################
## Section: calculate scaling factors using piRNAs or siRNAs and scale cpm using them and test if it works
######################################
######################################
stat.list = list.files(path = statDir, pattern = "*_cnt.typeHierarchy.txt", full.names = TRUE)
stats = NULL;
for(n in 1:length(stat.list))
{
  if(n ==1){
    stats = read.table(stat.list[n], sep = "\t", header=TRUE)
  }else{
    test = read.table(stat.list[n], sep = "\t", header=TRUE)
    stats = data.frame(stats, test[, -1])
  }
}
stats = stats[which(!is.na(stats$type)==TRUE), ]
rownames(stats) = stats[, 1]
stats = stats[, -1]

colnames(stats) = sapply(colnames(stats), function(x) gsub('count.', '', x), USE.NAMES = FALSE)

## need to merge again the techinical replicates for N2
if(Merge.techinical.replicates.N2){
  rep.technical = list(c("57751", "57753"), c("57752", "57754"))
  for(n in 1:length(rep.technical))
  {
    index = c()
    for(id in rep.technical[[n]])
    {
      #print(id)
      index = c(index, which(colnames(stats)==id))
    }
    
    #design$SampleID[index[1]] = paste0(design$SampleID[index], collapse = ".")
    ss = apply(stats[, index], 1, function(x) sum(x, na.rm = TRUE))
    stats[, (index[1])] = ss;
    colnames(stats)[index[1]] = paste0(colnames(stats)[index], collapse = ".")
      #paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  design$SampleID[index[1]])
    #design = design[-index[-1], ]
    stats = stats[, -index[-1]]
  }
} 

mm = match(design.matrix$SampleID, colnames(stats))

stats = stats[, mm]
colnames(stats) = paste0(design.matrix$genotype, "_", design.matrix$tissue.cell, "_", design.matrix$treatment, "_", design.matrix$SampleID)
stats = data.frame(t(stats))

source('RNAseq_Quality_Controls.R')
#pairs(stats, lower.panel=NULL, upper.panel=panel.fitting)

plot(library.sizes, stats$piRNA, log = 'xy')
plot(library.sizes, stats$siRNA, log = 'xy')

plot(stats$piRNA, stats$ncRNA, log='xy')
abline(0, 1, lwd=2.0, col='red')

#sizefactors.piRNA = stats$piRNA/median(stats$piRNA) 
#sizefactors.siRNA = stats$siRNA/median(stats$siRNA)
#sizefactors = (sizefactors.piRNA + sizefactors.siRNA) /2
sizefactors = stats$piRNA;
cpm.piRNA = countData
for(n in 1:ncol(cpm.piRNA))
{
  cpm.piRNA[,n] = countData[,n]/sizefactors[n]*10^6
}

## average the biological replicates
source("miRNAseq_functions.R")
cpm.piRNA.mean.rep = average.biological.replicates(cpm.piRNA)

## remove batch effect by scaling the untreated samples using N2 as the reference
cpm.piRNA.batch.corrected = remove.batch.using.N2.untreated(cpm.piRNA.mean.rep)

cpm.piRNA.batch.corrected[which(rownames(cpm.piRNA.batch.corrected)=='lsy-6'), grep('treated', colnames(cpm.piRNA.batch.corrected))]

## substract the background, namely the N2 in treated treatment
kk = grep('_treated', colnames(cpm.piRNA.batch.corrected))
index.N2 = intersect(grep('N2', colnames(cpm.piRNA.batch.corrected)), kk)

expression = cpm.piRNA.batch.corrected[, setdiff(kk, index.N2)] - cpm.piRNA.batch.corrected[, index.N2]   
expression[which(expression<0)] = 0

enriched.list = read.table(file = paste0(resDir, "/Enrichment_Matrix_13samples_66genes_with_clusters_for_neuronClasses.txt"), sep = "\t", 
                           header = TRUE, row.names = 1)
enriched.list = colnames(enriched.list)
enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)
mm = match((enriched.list), rownames(expression))

expression.sel = t(expression[mm, ])

library("pheatmap")
library("RColorBrewer")

pdfname = paste0(resDir, "/heatmap_ExpreMatrix_piRNAnormalization_for_12samples_66genes_with_clusters_for_neuronClasses", ".pdf")
pdf(pdfname, width=16, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

pheatmap(log2(expression.sel+1), cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=TRUE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

expression = data.frame(expression)
plot(expression$WT_Dopaminergic.neurons_treated, expression$WT_Ciliated.sensory.neurons_treated, log='xy')
abline(0, 1, lwd=2.0, col='red')

plot(expression$WT_mechanosensory.neurons_treated, expression$WT_unc.86.expressing.neurons_treated, log='xy')
abline(0, 1, lwd=2.0, col='red')

dev.off()

if(Save.Processed.Tables)
{
  write.table(expression.sel, file = paste0(resDir, "/Expression_Matrix_select_12samples_66genes_with_clusters_for_neuronClasses.txt"), 
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(expression, file = paste0(resDir, "/Expression_Matrix_select_12samples_allgenes_with_clusters_for_neuronClasses.txt"), 
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  write.table(cpm.piRNA.batch.corrected, file = paste0(resDir, "/Expression_Matrix_piRNA_normalization_average_replicates_remove_batch_12samples_allgenes_with_clusters_for_neuronClasses.txt"), 
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  
}




