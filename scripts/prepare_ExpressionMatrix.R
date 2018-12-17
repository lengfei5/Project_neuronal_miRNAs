##########################################################################
##########################################################################
## Project: Chiara's neuron class-specific miRNA expression  
## Script purpose: to test piRNA normalization for the expression matrix in the decomvolution analysis
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed May 16 13:58:52 2018
##########################################################################
##########################################################################
RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")
statDir = "../data/normalized_piRNAs"
version.table = "miRNAs_neurons_v1_2018_03_07"

resDir = "../results/tables_for_decomvolution"
tabDir = paste0(resDir, "/tables/")
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)
Save.Processed.Tables = TRUE
version.analysis = "neuronal_miRNAs_20181128"

calculate.sizeFactors.for.piRNAs = FALSE
Filter.lowly.expressed.using.predefined.miRNA.list = TRUE;

######################################
######################################
## Section: load count tables for miRNAs and prepare statistics for piRNA and siRNAs for normalization
######################################
######################################
load(file = paste0('../results/miRNAs_neurons_enrichment/Rdata/Design_Raw_readCounts_', version.table, '.Rdata'))
source("miRNAseq_functions.R")

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
    colnames(all)[(index[1]+1)] = paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  
                                         design$SampleID[index[1]])
    design = design[-index[-1], ]
    all = all[, -(index[-1]+1)]
  }
} 


####################
## calculate cpm and select only expressed miRNAs and start to use the gene names instead of arms 
####################
raw = as.matrix(all[, -1])
raw[which(is.na(raw))] = 0
library.sizes = apply(raw, 2, sum)
raw = floor(raw)
rownames(raw) = all$gene

if(Filter.lowly.expressed.using.predefined.miRNA.list){
  ## filter mir star and also the genes with low counts
  source("miRNAseq_functions.R")
  expressed.miRNAs =  find.expressed.mature.miRNA.using.cpm.threshold(raw[, which(design$treatment == "untreated")], 
                                                                      cpm.threshold = 0.1)
  expressed.miRNAs = data.frame(expressed.miRNAs[, c(1:3)], stringsAsFactors = FALSE)
  
}else{
  # filter lowly expressed miRNA with list of predefined miRNAs that were identified using all untreated samples 
  list.expressed = read.csv(paste0("../data/list_expressed_miRNAs_using_Untreated_samples_Henn1_mutant_WT_all_cpm_10.csv"), header = TRUE, as.is = c(1, 2))
  # prepare old llist
  list.expressed = find.mature.ones.for.prefixed.expressed.miRNAs(list.expressed)
  expressed.miRNAs = data.frame(list.expressed[, c(1:3)], stringsAsFactors = FALSE)
  
}
expressed.miRNAs = expressed.miRNAs[expressed.miRNAs$mature,]
mm = match(expressed.miRNAs$miRNA, all$gene)
raw = raw[mm, ]
rownames(raw) = expressed.miRNAs$gene

###############################
# filter the samples unrelevant (non-neuron samples, henn1-mutant background)
###############################
kk =  unique(c(which(design$tissue.cell=="Glial.cells" | design$tissue.cell == "CEPsh" | 
             design$genotype == "henn1.mutant"), grep("L3", design$tissue.cell)))

if(length(kk)>0){
  design.matrix = design[-kk, ]
  countData = raw[, -kk]
  library.sizes = library.sizes[-kk]
}

######################################
######################################
## Section: piRNA normalization 
# 1) import and reform the piRNAs, siRNAs statistics for all samples
# 2) try to calculate the sizefactors from individual piRNA 
######################################
######################################
load(file = "../results/tables_for_decomvolution/Rdata/rawCounts_table_for_piRNAs.Rdata") # piRNA counts
rownames(piRNA.counts) = piRNA.counts$gene
piRNA.counts = as.matrix(piRNA.counts[, -1])
piRNA.counts[which(is.na(piRNA.counts))] = 0

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
  source("miRNAseq_functions.R")
  
  rep.technical = list(c("57751", "57753"), c("57752", "57754"))
  stats = Merge.techinical.replicates(stats = stats, rep.technical = rep.technical)
  piRNA.counts = Merge.techinical.replicates(stats = piRNA.counts, rep.technical = rep.technical)
}

mm = match(design.matrix$SampleID, colnames(stats))
stats = stats[, mm]
colnames(stats) = paste0(design.matrix$genotype, "_", design.matrix$tissue.cell, "_", design.matrix$treatment, "_", design.matrix$SampleID)
stats = data.frame(t(stats))

mm = match(design.matrix$SampleID, colnames(piRNA.counts))
piRNA.counts = piRNA.counts[, mm]
colnames(piRNA.counts) = paste0(design.matrix$genotype, "_", design.matrix$tissue.cell, "_", design.matrix$treatment, "_", design.matrix$SampleID)

save(stats, piRNA.counts, countData, design.matrix, library.sizes, 
     file = paste0(RdataDir, 'neuronalClasses_samples_countTables_piRAN_siRNA_stats_', version.table, '.Rdata'))

########################################################
########################################################
# Section: Normalization with piRNAs
# double check sample quality
# calculate scaling factors using piRNAs or siRNAs and scale cpm using them 
# check again the sample qualities after normalization 
########################################################
########################################################
Compare.piRNA.siRNA.spikeIns.as.scaling.factors = FALSE
if(Compare.piRNA.siRNA.spikeIns.as.scaling.factors){
  
  load(file = paste0(RdataDir, 'neuronalClasses_samples_countTables_piRAN_siRNA_stats_', version.table, '.Rdata'))
  
  pdfname = paste0(resDir, "/Compare_piRNAs_siRNAs_spikeIns_for_scalingFactorinNormalization", version.analysis,  ".pdf")
  pdf(pdfname, width=20, height = 14)
  par(cex =0.7, mar = c(5,5,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  # common quality control normalized by DESeq2 
  source("RNAseq_Quality_Controls.R")
  Check.RNAseq.Quality(read.count=countData, design.matrix = design.matrix[, c(1,3,5)], lowlyExpressed.readCount.threshold = 0)
  
  source("miRNAseq_functions.R")
  Compare.piRNA.siRNA.spikeIns.for.scaling.factors(library.sizes, stats, countData, design.matrix)
  
  #Compare.piRNA.librarysize.vs.sizeFactors.vs.spikeIns(library.sizes, stats, piRNA.counts, countData, design.matrix)
  
  dev.off()
  
  pdfname = paste0(resDir, "/QC_for_piRNAs_countTable_", version.analysis,  ".pdf")
  pdf(pdfname, width=20, height = 14)
  par(cex =0.7, mar = c(5,5,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  Check.RNAseq.Quality(read.count=piRNA.counts, design.matrix = design.matrix[, c(1, 3,5)], lowlyExpressed.readCount.threshold = 20)
  
  dev.off()
  
  pdfname = paste0(resDir, "/Compare_piRNAs_librarySize_vs_sizeFactors_spikeIns_", version.analysis,  ".pdf")
  pdf(pdfname, width=20, height = 14)
  par(cex =0.7, mar = c(5,5,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  source("miRNAseq_functions.R")
  piRNA.sizeFctors = calculate.sizeFactors4piRNAs(read.count = piRNA.counts, design.matrix = design.matrix[, c(1, 3, 5)], 
                                                  lowlyExpressed.readCount.threshold = 100)
  
  plot(piRNA.sizeFctors, stats$piRNA, log = 'xy', xlab = "size factors", ylab = 'piRNA library size', cex = 1.5)
  Compare.piRNA.siRNA.spikeIns.for.scaling.factors(library.sizes, stats, countData, design.matrix, 
                                                                      piRNA.sizeFctors = piRNA.sizeFctors)
  
  dev.off()

}

####################
## Normalize the data using piRNAs
####################
#sizefactors.piRNA = stats$piRNA/median(stats$piRNA) 
#sizefactors.siRNA = stats$siRNA/median(stats$siRNA)
#sizefactors = (sizefactors.piRNA + sizefactors.siRNA) /2
sizefactors = stats$piRNA;
cpm.piRNA = countData
for(n in 1:ncol(cpm.piRNA))
{
  cpm.piRNA[,n] = countData[,n]/sizefactors[n]*10^6
}

save(stats, countData, design.matrix, cpm.piRNA, file = paste0(RdataDir, 'piRAN_siRNA_stats_counTables_cpm.piRNA_', version.table, '.Rdata'))

######################################
######################################
## Section: remove batch effects and check the expression matrix
######################################
######################################
load(file = paste0(RdataDir, 'piRAN_siRNA_stats_counTables_cpm.piRNA_', version.table, '.Rdata'))

Use.ComBat.for.batch.correction = TRUE
## to test if one rab-3 replicates should be removed
Remove.pan.neurons.samples.71822.71823.as.outliers = TRUE

###############################
# choose the method of batch correction 
###############################
source("miRNAseq_functions.R")
design.matrix$batch = c(rep(1, 4), rep(2, 2), rep(c(3:14), each=4), rep(17, 4), rep(15, 2), rep(16, 2), rep(18, 8))

if(Remove.pan.neurons.samples.71822.71823.as.outliers)
{
  kk = match(c(71822:71823), design.matrix$SampleID)
  #cpm.piRNA.bc.meanrep = average.biological.replicates(cpm.piRNA.bc[, -kk])
  design.matrix = design.matrix[-kk, ]
  cpm.piRNA = cpm.piRNA[, -kk]
}

#method.sel = 'linear.model'
## remove batch effect by scaling the untreated samples using N2 as the reference
pdfname = paste0(resDir, "/Check_batch_difference_",  version.analysis, ".pdf")
pdf(pdfname, width=12, height = 12)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
cpm.piRNA.bc.my = remove.batch.using.N2.untreated(cpm.piRNA, design.matrix, method = 'linear.model')
dev.off()

cpm.piRNA.bc.limma = remove.batch.using.N2.untreated(cpm.piRNA, design.matrix, method = 'limma')
cpm.piRNA.bc.combat = remove.batch.using.N2.untreated(cpm.piRNA, design.matrix, method = 'combat')

#cpm.piRNA.bc.ratio = remove.batch.by.logratios(cpm.piRNA, design.matrix)

cpm.piRNA.bc = cpm.piRNA.bc.combat

#if(Use.ComBat.for.batch.correction){
#}else{ cpm.piRNA.bc = cpm.piRNA.bc.my }

## double check the batch correction
source("RNAseq_Quality_Controls.R")
pdfname = paste0(resDir, "/Check_batchRemoval_for_rab-3_rgef-1",  version.analysis, ".pdf")
pdf(pdfname, width=12, height = 12)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

jj = grep("WT_Pan.neurons_untreated_", colnames(cpm.piRNA.bc))
kk = grep("WT_Pan.neurons_treated_", colnames(cpm.piRNA.bc))

par(mfrow=c(1, 1))
plot.pair.comparison.plot(cpm.piRNA.bc[, jj], main = paste0("untreate pan-neurons" ,"-- pairwise comparison for piRNA-normalization"))
plot.pair.comparison.plot(cpm.piRNA.bc[, kk], main = paste0("treate pan-neurons" ,"-- pairwise comparison for piRNA-normalization"))

dev.off()

pdfname = paste0(resDir, "/Check_piRNA_normalization_batchRemoval_",  version.analysis, ".pdf")
pdf(pdfname, width=18, height = 6)
par(cex =0.7, mar = c(6,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
#par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))
source("miRNAseq_functions.R")

Test.piRNA.normalization.batch.removal(cpm.piRNA, design.matrix)
#Test.piRNA.normalization.batch.removal(cpm.piRNA.bc.my, design.matrix)
Test.piRNA.normalization.batch.removal(cpm.piRNA.bc.limma, design.matrix)
Test.piRNA.normalization.batch.removal(cpm.piRNA.bc.combat, design.matrix)

dev.off()

## average the biological replicates
source("miRNAseq_functions.R")
cpm.piRNA.bc.meanrep = average.biological.replicates(cpm.piRNA.bc)
cpm.piRNA.bc.meanrep.log2 = average.biological.replicates(log2(cpm.piRNA.bc))

save(cpm.piRNA.bc, 
     cpm.piRNA.bc.meanrep, 
     cpm.piRNA.bc.meanrep.log2,
     design.matrix, 
     file = paste0(RdataDir, 'piRANormalized_cpm.piRNA_batchCorrectedCombat_reAveraged_', version.table, '.Rdata'))


########################################################
########################################################
# Section : calibrate promoter methylation efficiencies
# after piRNA normalization and bacth correction using untreated samples
# The background for treated samples of differetn promoters are different (the assumption)
# we are using the non-enriched miRNAs to correct this bias, which are resulted from the background composition, promoter mythelation efficiencies.
# hopefully this will not change too much our initial data and the background correction will be roughly proprotional to the sample sizes
########################################################
########################################################
load(file = paste0(RdataDir, 'piRANormalized_cpm.piRNA_batchCorrectedCombat_reAveraged_', version.table, '.Rdata'))
####################
## Here we decided to use the piRNA normalization and correct the batch using ComBat 
####################
Try.Calibrate.promoter.methylation.efficiency = TRUE
if(Try.Calibrate.promoter.methylation.efficiency){
  pdfname = paste0(resDir, "/Select_panNeurons_BEFORE_background_calibration_Samples_check_examples_",  version.analysis, ".pdf")
  pdf(pdfname, width=12, height = 10)
  par(cex =0.7, mar = c(6,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  source("miRNAseq_functions.R")
  Compare.pan.neuron.vs.other.five.samples.And.check.miRNA.examples(cpm.piRNA.bc, design.matrix)
  
  dev.off()
  
  pdfname = paste0(resDir, "/Select_panNeurons_AFTER_background_calibration_Samples_check_examples_",  version.analysis, ".pdf")
  pdf(pdfname, width=12, height = 10)
  par(cex =0.7, mar = c(6,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  source("miRNAseq_functions.R")
  xx = calibrate.promoter.methylation.efficiency(cpm.piRNA.bc, design.matrix)
  cpm.piRNA.bc.prot = 2^xx
  
  Compare.pan.neuron.vs.other.five.samples.And.check.miRNA.examples(cpm.piRNA.bc.prot, design.matrix)
  #(cpm.piRNA.bc.meanrep.log2[mm, grep("_treated", colnames(cpm.piRNA.bc.meanrep.log2))])
  #Compare.pan.neuron.vs.other.five.samples(cpm.piRNA.bc)
  dev.off()
  
}

save(cpm.piRNA.bc.prot,
     cpm.piRNA.bc, 
     cpm.piRNA.bc.meanrep, 
     cpm.piRNA.bc.meanrep.log2,
     design.matrix, 
     file = paste0(RdataDir, 'piRANormalized_cpm.piRNAnorm_batchCorrectedCombat_calibratedProtEff_', version.table, '.Rdata'))

