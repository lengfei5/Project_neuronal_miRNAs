##########################################################################
##########################################################################
## Project: Chiara's neuron class-specific miRNA expression
## Script purpose: this script is the main script for the deconvolution anlaysis which takes the fraction matrix, expression matrix prepared by 
## prepare_FractionMatrix.R and prepare_ExpressionMatrix.R 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Fri May 25 08:00:23 2018
##########################################################################
##########################################################################
library("pheatmap")
library("RColorBrewer")
source('miRNAseq_functions.R')
RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")

version.ExprsMatrix = "miRNAs_neurons_v1_2018_03_07"
version.Fraction.Matrix = "_miRNAs_neurons_20180525"
version.EnrichscoreMatrix = "20180506"

version.analysis = "_20180920"

resDir = "../results/decomvolution_results/"
if(!dir.exists(resDir)) dir.create(resDir)
#testDir = paste0(resDir, "deconv_test_09_20_gcdnet_tuning_global_lambda2_cv/")
#if(!dir.exists(testDir)) dir.create(testDir)

Data.complete = TRUE
fitting.space = "log2" ## linear or log2 transformed for expression matrix

add.background.sample.in.fitting.linear.space = TRUE
Use.mergedFractionMatrix = TRUE 
Use.mergedExpressionMatrix = FALSE # group the genes if they show similar gene expression pattern
Manually.unifiy.sample.names.forMatrix = TRUE
Check.ProprotionMatrix.ExpressionMatrix = FALSE

######################################
######################################
## Section: load fraction matrix A and merge similar neuron classes
######################################
######################################
if(Use.mergedFractionMatrix){
  load(file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix_plus_mergedFractionMatrix", 
                     version.Fraction.Matrix, ".Rdata"))
  proportions = newprop;
  
}else{
  load(file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix", 
                     version.Fraction.Matrix, ".Rdata"))
  source('miRNAseq_functions.R')
  
  pdfname = paste0(resDir, "/heatmap_for_merging_proportionaMatrix.pdf")
  pdf(pdfname, width=15, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  
  newprop = proportions.matrix.merging.neuronClass(proportions)
  
  dev.off()
  
  save(newprop, mapping, neurons, proportions, 
       file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix_plus_mergedFractionMatrix", 
                     version.Fraction.Matrix, ".Rdata"))
}

######################################
######################################
## Section: enrichment score matrix and expression matrix
######################################
######################################
if(!Use.mergedExpressionMatrix){
  load(file = paste0(RdataDir, 'piRANormalized_cpm.piRNA_batchCorrectedCombat_reAveraged_', version.ExprsMatrix, '.Rdata'))
  jj = grep('_untreated', colnames(cpm.piRNA.bc.meanrep))
  total = apply(cpm.piRNA.bc.meanrep[, jj], 1, median)
  xx = data.frame(total, cpm.piRNA.bc.meanrep[, -jj])
  ncs = sapply(colnames(xx)[-c(1:2)], function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
  ncs = sapply(ncs, function(x) gsub("*.neurons", "", x), USE.NAMES = FALSE)
  colnames(xx) = c('whole.body', 'background', ncs) 
  
  ####################
  ## here we transform the gene expression by e'= (expression-background)/background 
  # then e' = 0 if e'<0 or e'<1; 
  # this transformation is to scale the range for each gene so that all data fall into the same range and e' should still follow the positive constrain 
  # meanwhile using ratio between expression and background to filter non-expressed ones  
  ####################
  expression = xx[, -c(1)] # ignore the gene expression in the whole body
  
  
}else{
  source('miRNAseq_functions.R')
  newExprM = expressionMatrix.grouping(xx) 
}

########################################################
########################################################
# Section: select the miRNAs to analyze 
# Match the samples for fraction matrix and expression matrix
# especailly for the fitting in the linear scale
########################################################
########################################################
###############################
# select samples to use and also match sample orders in matrix A and Y 
# mannually unify the names of each sample A and Y 
###############################
if(Manually.unifiy.sample.names.forMatrix){
  rownames(proportions) = c("Dopaminergic","Serotonergic","Glutamatergic", "Cholinergic", "GABAergic",  "Ciliatedsensory",
                            "Mechanosensory", "unc.86", "Pharyngeal", "ceh.14", "unc.42", "unc.3","ASE", "Pan.neurons")
  expression = t(expression)
  rownames(expression) = c( "background", "ASE", "Serotonergic", "Dopaminergic", "Glutamatergic", "Ciliatedsensory",
                            "GABAergic", "Mechanosensory", "unc.86", "Pharyngeal", "Cholinergic", "ceh.14", "unc.42", 
                            "Pan.neurons", "unc.3")
  
  ###############################
  # if fitting in linear scale, then the background sample is used and the background is added in proportion matrix as well
  # if fitting in log2 scale, then the log2(sample/background) is used; in this case, the coefficients must be also positive
  ###############################
  if(fitting.space == "linear"){
    if(add.background.sample.in.fitting.linear.space){
      yy = proportions;
      yy = rbind(rep(0, ncol(proportions)), proportions)
      yy = cbind(rep(1, nrow(yy)), yy)
      colnames(yy) = c("C0", colnames(proportions))
      rownames(yy) = c("background", rownames(proportions))
      proportions = yy
    }
  }else{
    
    yy = expression;
    # log2 (expresion/background)
    for(n in 1:nrow(yy)) yy[n, ] = log2(as.numeric(yy[n, ])/as.numeric(expression[1,]))
    
    yy = yy[-1, ]
    expression = yy;
    
  }
  
  save(expression, proportions, file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                                              version.Fraction.Matrix,"_", version.ExprsMatrix, "fitting_scale", 
                                              fitting.space,".Rdata"))
  #index.sel = c(13, 2, 1, 3, 6, 5, 7, 8, 9, 4, 10, 11)
  #proportions.sel = proportions[index.sel, ]
}

###############################
# select genes of intereset and match sample orders
###############################
load(file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                     version.Fraction.Matrix,"_", version.ExprsMatrix, "fitting_scale", 
                     fitting.space,".Rdata"))

load(file = paste0(RdataDir, "Enrichscores_Matrix_13samples_selected_and_all_genes_", version.EnrichscoreMatrix, ".Rdata"))
enriched.list = colnames(enrich.matrix.sel)
#enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)

mm = match((enriched.list), colnames(expression))
expression.sel = expression[, mm]
#expression.sel = log2(expression.sel)

sels = match(rownames(expression.sel), rownames(proportions))
proportions.sel = proportions[sels, ]

if(Data.complete){
  jj = match(rownames(proportions), rownames(proportions.sel))
  proportions.sel = proportions.sel[jj, ]
  expression.sel = expression.sel[jj, ]
}

####################
## double check the proprotion matrix and expression matrix
## match the sample order in the proprotion matrix and expression matrix 
## now manually (to change)
####################
if(Check.ProprotionMatrix.ExpressionMatrix){
  xx = proportions.sel;
  xx[which(xx>0)] = 1
  
  yy = expression.sel;
  #yy = expression.sel[c(3, 2, 4, 10, 6, 5, 7, 8, 9, 11, 12, 1), ];
  #yy = proportions[c(index.sel, 12, 14), ]
  
  library("pheatmap")
  library("RColorBrewer")
  
  pdfname = paste0(resDir, "/Heatmap_Proportiona_Expression_Matrix", version.analysis, ".pdf")
  pdf(pdfname, width=15, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  # par(mfcol=c(1, 1))
  
  if(fitting.space == 'linear') {
    logaxis = 'xy';
    yy = t(log2(t(yy)/yy[which(rownames(yy)=='background'), ]));
    
    xx = xx[-1, -1];
    yy = yy[-1,  ]
    
  }else{logaxis = ''}
  
 
  pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=TRUE, 
           color = c("lightgray", "blue"), legend = FALSE)
  
  pheatmap(yy, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, 
           cluster_cols=TRUE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))
  
  ## double check the expression matrix
  par(mfrow = c(1, 1))
  mm = match(c("Cholinergic", "Glutamatergic",  "GABAergic",  "Dopaminergic", "Serotonergic"), rownames(yy))
  plot(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum), xlab = "Pan.neurons", ylab = "sum of Cho, Glut, GABA, Dop and Ser")
  abline(0, 1, col="red", lwd=2.0)
  #abline(h=1, col="darkgray", lwd=2.0)
  
  text(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum), labels = colnames(yy), cex = 0.8,
       pos = 1, offset = 0.4)
  
  par(mfrow = c(1, 2))
  plot(t(expression.sel[match(c("Dopaminergic", "Ciliatedsensory"), rownames(expression.sel)), ]), log=logaxis)
  abline(0, 1, lwd=2.0, col='red')
  plot(t(expression.sel[match(c("Mechanosensory",  "unc.86"), rownames(expression.sel)), ]), log=logaxis)
  abline(0, 1, lwd=2.0, col='red')
  
  dev.off()
  
}

######################################
######################################
# Section: process the table and run the glmnet
# the glmnet will be run for each gene, because the group lasso is not desirable
######################################
######################################
require(glmnet)
require(gcdnet)
library("pheatmap")
library("RColorBrewer")

x=as.matrix(proportions.sel)
y = as.matrix(expression.sel)

x = x >0 

## force value to be zero if they are lower than the background
y[y<0] = 0

Example2test = c("lsy-6", "mir-791", "mir-790", "mir-793",  "mir-792","mir-1821", "mir-83", "mir-124")
jj2test = match(Example2test, colnames(y))
y = y[, jj2test[which(!is.na(jj2test)==TRUE)]]

res = matrix(NA, nrow = ncol(x), ncol = ncol(y)) 
colnames(res) = colnames(y)
rownames(res) = colnames(x)
#x.ms = apply(x, 2, sum)
#x = x[, which(x.ms>0)]
#x = x>0


########################################################
########################################################
# Section: gcdnet with global lambda2 parameter
########################################################
########################################################
TEST.gcdnet.global.lambda2 = TRUE
if(TEST.gcdnet.global.lambda2){
  
  testDir = paste0(resDir, "decon_TEST/deconv_test_09_21_gcdnet_tuning_global_lambda2")
  if(!dir.exists(testDir)) dir.create(testDir)
  
  require(gcdnet)
  source("select_tuningParams_elasticNet.R")
    
  cols = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100)
  
  lambda2s = c(seq(0.01, 0.09, by=0.01), seq(0.1, 0.5, by=0.05))
  lambda = 10^seq(-4, 2, length.out = 200)
  
  for(method in c("cv.lambda.1se"))
  {
    pdfname = paste0(testDir, "deconv_test", version.analysis, 
                     "_fitting.", fitting.space, "_global_lambda2_",lambda2, "_method4crit_", method, ".pdf")
    
    pdf(pdfname, width=16, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    par(mfrow=c(1, 1))
    
    for(lambda2 in lambda2s)
    {
      # lambda2 = 0.2
      cat("lambda2 --", lambda2, "----------\n")
      #alpha = 0.7 # 0.001 is the default valu
      for(n in 1:ncol(y))
      {
        # n = 1; lambda2 = lambda2s[1];
        source("select_tuningParams_elasticNet.R")
        main = paste0(" (lambda2=", signif(lambda2, d=2), ")::", colnames(y)[n]);
        myCoefs = select.tuning.parameters.for.gcdnet(x, y[,n], lambda = lambda, lambda2=lambda2, nfold =10, method = method, omit.BIC = TRUE);
        res[,n] = as.numeric(myCoefs$coefficients)
        
        #par(mfrow= c(1,1))
        #plot(cv.fit, main = colnames(y)[n])
        
        #cat("---", colnames(res)[n], "---\n", 
        #    paste0(rownames(res)[which(res[,n]>0)], collapse = "\n"), "\n")
        
        #coef(cv.fit, s = "lambda.min")
        #myCoefs[which(myCoefs != 0 ) ]               #coefficients: intercept included
        #myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ] #feature names: intercept included
      }
      pheatmap((res), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
               cluster_cols=FALSE, main = paste0("lambda2 = ", lambda2, " - method : ", method), na_col = "gray30",
               color = cols)
      
    }
  
    dev.off() 
  }
   
}

########################################################
########################################################
# Section: glmnet with global alpha parameter or gene-specific alpha parameters
########################################################
########################################################
require(glmnet)

TEST.glmnet.gene.specific.alpha = TRUE

Methods2test = c("cv.lambda.1se")

alphas = c(0.01, 0.02, 0.05,seq(0.1, 1, by= 0.05))
#alphas = c(0.1)
lambda = 10^seq(-4, 2, length.out = 200)

source("select_tuningParams_elasticNet.R")

if(TEST.glmnet.gene.specific.alpha) {testDir = paste0(resDir, "decon_TEST/deconv_test_09_21_glmnet_tuning_gene_specific_alpha");
}else{testDir = paste0(resDir, "decon_TEST/deconv_test_09_21_glmnet_tuning_global_alpha"); }
if(!dir.exists(testDir)) dir.create(testDir)

for(method in Methods2test)
{
  pdfname = paste0(testDir, "/deconv_test", version.analysis, 
                   "_fitting.", fitting.space, 
                   "_glmnet_global_alpha_method_select_tuning_parameters_", method,
                   ".pdf")
  
  pdf(pdfname, width=16, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  
  run.glmnet.select.tuning.parameters(x, y, alphas = alphas, method = method, Gene.Specific.Alpha = TEST.glmnet.gene.specific.alpha);
  
  dev.off()
  
}


