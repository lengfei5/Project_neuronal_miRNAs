##################################################
##################################################
## Project: Chiara's neuroclass specific miRNA expression
## Script purpose: to make the fraction matrix A in the formular AX = Y 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Thu Apr 26 15:40:49 2018
##################################################
##################################################
find.neuron.class = function(x) ## NOT used here
{
  library(stringr)
  if(length(grep("^DD", x))>0) {
    x = 'DD'
  }else{
    if(length(grep("^DA", x))>0){
      x = "DA"
    }else{
      if(length(grep("^DB", x))>0){
        x = "DB"
      }else{
        y = gsub('DR$', '', as.character(x))
        if(str_length(y)>1) x = y
        y = gsub('DL$', '', as.character(x))
        if(str_length(y)>1) x = y
        y = gsub('VL$', '', as.character(x))
        if(str_length(y)>1) x = y
        y = gsub('VR$', '', as.character(x))
        if(str_length(y)>1) x = y
        x = gsub('L$', '', as.character(x))
        x = gsub('R$', '', as.character(x))
      }
    }
  }
  return(x)
}

Save.Processed.Tables = FALSE
version.Fraction.Matrix = "_miRNAs_neurons_20180525"

library(openxlsx)
resDir = "../results/tables_for_decomvolution"
RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

####################
## save the sample information and promoters
####################
aa = read.xlsx("../annot_neurons/Matrix4deconvolution_Jingkui_April_25th_L1.xlsx", sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, 
               skipEmptyCols = TRUE)

### clean the matrix 
bb = aa[, -c(2, 3, 4, 5, 18, 19)]
samples = colnames(bb)[-1]
samples = c(samples, 'ASE', "Pan-neurons")
promoters = c('dat-1', 'tph-1', 'eat-4', 'unc-17', 'unc-47', 'osm-5', 'mec-3', 'unc-86', 'ehs-1', 'ceh-14', 'unc-42', 'unc-3', 'ASE', 'rab-3')
neurons = sapply(samples, function(x) unlist(strsplit(as.character(x), "[.]"))[1])

mapping = data.frame(samples, promoters, neurons, stringsAsFactors = FALSE)
rownames(mapping) = c(1:nrow(mapping))
rownames(mapping) = c(1:nrow(mapping))
colnames(mapping) = c('sampleInfo', 'promoters', 'neuron.sampling')
#colnames(bb)[-1] = mapping$promoters
if(Save.Processed.Tables){
  write.table(mapping, file = paste0(resDir, "/sample_infos_promoters.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

####################
## Construct proportion matrix  for neuron classes  
####################
classes = read.xlsx("../annot_neurons/Neurons_vs_neuronal_classes.xlsx", sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, skipEmptyCols = TRUE)
aa = read.xlsx("../annot_neurons/Matrix4deconvolution_Jingkui_April_25th_L1.xlsx", sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, 
               skipEmptyCols = TRUE)
bb = aa[, -c(2, 3, 4, 5, 18, 19)]
neurons = bb;

## here add ASE and Pan-neurons 
neurons$ASE = NA; 
neurons$ASE[grep('ASE', neurons$Neurons.in.L1)] = "ASE"
neurons$Pan.neurons = "pan.neurons";
neurons$Pan.neurons[grep('CAN', neurons$Neurons.in.L1)] = NA

neurons = neurons[, c(2:15, 1)]
mm = match(neurons$Neurons.in.L1, classes$`222.NEURONS.IN.L1`)
neurons$neuronClass.in.L1 = classes$NEURONAL.CLASS.FOR.EACH.NEURON[mm]
jj = match(colnames(neurons)[c(1:14)], mapping$sampleInfo)
colnames(neurons)[c(1:14)] = mapping$neuron.sampling[c(1:14)]
colnames(neurons)[15:16] = c('neurons.L1', 'neuronClass.L1')

presence = as.matrix(neurons[, c(1:14)])
presence[which(!is.na(presence))] = 1
presence[which(is.na(presence))] = 0
presence = data.matrix(presence)

neurons = data.frame(presence, neurons[, -c(1:14)])
write.table(neurons, file = paste0(resDir, "/Presence_Absence_neurons_for_each_sample.txt"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

neurons = read.table(file = paste0(resDir, "/Presence_Absence_neurons_for_each_sample.txt"), sep = "\t", header = TRUE, as.is = c(15, 16))

proportions = matrix(NA, ncol = 14, nrow = length(unique(neurons$neuronClass.L1)))
colnames(proportions) = colnames(neurons)[c(1:14)];
rownames(proportions) = unique(neurons$neuronClass.L1)

for( n in 1:nrow(proportions))
{
  # n = 1;
  cat(n, '--', rownames(proportions)[n], "\n")
  jj = which(neurons$neuronClass.L1 == rownames(proportions)[n])
  if(length(jj)>1){
    test = apply(((neurons[jj, c(1:14)])), 2, sum);
    #test[which(test>0)] = 1
  }else{
    test = neurons[jj, c(1:14)]
  }
  proportions[n, ] = as.numeric(test);
}

proportions = t(proportions)
ss = apply(proportions, 2, sum)

if(Save.Processed.Tables)
{
  save(mapping, neurons, proportions, 
       file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix", version.Fraction.Matrix, ".Rdata"))
  
  write.table(proportions, file = paste0(resDir, "/Proportions_Matrix_14_samples_96_neuronClasses.txt"), sep = "\t", col.names = TRUE, row.names = TRUE, 
              quote = FALSE)
}

##################################################
##################################################
## Section: Check the proportion matrix and check correlation between eahc feature
##################################################
##################################################
library("pheatmap")
library("RColorBrewer")

pdfname = paste0(resDir, "/heatmap_for_proportion_matrix_14_samples_96_neuronClasses_with_clusters_for_neuronClasses", ".pdf")
pdf(pdfname, width=16, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

pheatmap(proportions, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=TRUE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

dev.off()


#nns = sapply(neurons, find.neuron.class)
#nns.uniq = unique(nns)


