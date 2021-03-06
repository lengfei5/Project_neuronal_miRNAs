##################################################
##################################################
## Project: Chiara's neuron-specific miRNA expression 
## Script purpose: to prepare the expression matrix in the deconvolution analysis and test the normalization using pi-RNAs
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue May 15 17:50:26 2018
##################################################
##################################################

dataDir = "../results/miRNAs_neurons_enrichment/tables"
RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")

version.analysis = "20181203"
Save.Processed.Tables = FALSE

resDir = "../results/tables_for_decomvolution"
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

##################################################
##################################################
## Section: construct the matrix using the enrichment analysis
# select the genes of interest: the ones showing enrichment in any of those samples
# the cutoff should be considered and specified
##################################################
##################################################
enrich.files = list.files(path = dataDir, pattern = paste0("*_Mature_miRNAs_neurons_v1_2018_03_07.csv"), full.names = TRUE, 
                          recursive = TRUE, include.dirs = TRUE)
#jj = grep("ASE.neurons_henn1.mutant", enrich.files)
#enrich.files = enrich.files[-jj]

nsamples = unlist(sapply(enrich.files, function(x) unlist(strsplit(basename(x), "_"))[4]))
genotypes = unlist(sapply(enrich.files, function(x) unlist(strsplit(basename(x), "_"))[5]))
promoters = unlist(sapply(enrich.files, function(x) unlist(strsplit(basename(x), "_"))[6]))

####################
## Construct the matrix using all samples (three promoters for pan-neurons (all henn1-mutant) and one ASE (WT) 
####################
enrich.matrix = NULL
for(n in 1:length(enrich.files))
{
  test = read.csv(enrich.files[n], header = TRUE, row.names = 1)
  if(n==1){
    enrich.matrix = data.frame(test$log2FoldChange, test$pvalue)
    rownames(enrich.matrix) = rownames(test)
  }else{
    mm = match(rownames(enrich.matrix), rownames(test))
    enrich.matrix = data.frame(enrich.matrix, test$log2FoldChange[mm], test$pvalue[mm])
  }
  colnames(enrich.matrix)[c((2*n-1),2*n)] = paste0(nsamples[n], "_", genotypes[n], "_",  promoters[n], "_", c('log2FC', 'pvalue'))
}

## determine candidates of interest 
## Log2 FC > 1 and p val< 0.01 for WT back; Log2 FC > 1 and p val< 0.001 for mutant
fc.cutoff = 1; 
pval.muant = 0.001;
pval.wt = 0.01;

candidates = rep(0, nrow(enrich.matrix))
names(candidates) = rownames(enrich.matrix)
for(n in 1:length(candidates))
{
  nb.enrich = 0
  for(m in 1:length(nsamples))
  {
    if(nsamples[m] == "Pan.neurons"){pval.cutoff = pval.muant;
    }else{pval.cutoff = pval.wt}
    
    if(enrich.matrix[n, (2*m-1)] >= fc.cutoff & enrich.matrix[n, 2*m]<pval.cutoff) nb.enrich = nb.enrich +1
  }
  
  cat(names(candidates)[n], "--", nb.enrich, '--fc.cutoff ', fc.cutoff, "--pval.cutoff ", pval.cutoff,  "\n")
  
  if(nb.enrich>0) candidates[n] = 1
}

##
## concatenate Pan-neurons
pan.ns = as.matrix(enrich.matrix[, grep("Pan.neurons", colnames(enrich.matrix))])
pan.ns.cat = data.frame(apply(pan.ns[, grep('_log2FC', colnames(pan.ns))], 1, mean), 
                        apply(pan.ns[, grep('_pvalue', colnames(pan.ns))], 1, min))
colnames(pan.ns.cat) = paste0('Pan.neurons_rab.3.rgef.1.unc.31', c('_log2FC', "_pvalue"))
enrich.matrix = data.frame(pan.ns.cat, enrich.matrix[, grep("Pan.neurons", colnames(enrich.matrix), invert = TRUE)])

enrich.matrix.sel = t(enrich.matrix[which(candidates>0), grep('_log2FC', colnames(enrich.matrix))])

####################
## save the enrichment matrix 
####################
library("pheatmap")
library("RColorBrewer")

pdfname = paste0(resDir, "/heatmap_for_13samples_66genes_with_clusters_for_neuronClasses", version.analysis, ".pdf")
pdf(pdfname, width=16, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

pheatmap(enrich.matrix.sel, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=TRUE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

dev.off()

if(Save.Processed.Tables)
{
  save(enrich.matrix.sel, enrich.matrix, 
       file = paste0(RdataDir, "Enrichscores_Matrix_13samples_selected_and_all_genes_", version.analysis, ".Rdata"))
  write.table(enrich.matrix.sel, file = paste0(resDir, "/Enrichment_Matrix_13samples_66genes_with_clusters_for_neuronClasses.txt"), 
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(enrich.matrix, file = paste0(resDir, "/Enrichment_Matrix_13samples_allgenes_with_clusters_for_neuronClasses_", version.analysis, ".txt"), 
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
}




