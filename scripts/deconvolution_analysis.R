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
RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")

version.ExprsMatrix = "miRNAs_neurons_v1_2018_03_07"
version.Fraction.Matrix = "_miRNAs_neurons_20180525"
version.EnrichscoreMatrix = "20180506"

version.analysis = "_20180802"

resDir = "../results/decomvolution_results"
if(!dir.exists(resDir)) dir.create(resDir)

fitting.space = "linear" ## linear or log2 transformed for expression matrix
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
  if(fitting.space == "linear")
  {
    expression = xx[, -c(1)] # ignore the gene expression in the whole body
  }else{
    expression = xx[, -c(1, 2)]
    expression = log2(expression/xx$background)
  }
    
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
  if(fitting.space == "linear"){
    expression = t(expression)
    rownames(expression) = c( "background", "ASE", "Serotonergic", "Dopaminergic", "Glutamatergic", "Ciliatedsensory",
                                  "GABAergic", "Mechanosensory", "unc.86", "Pharyngeal", "Cholinergic", "ceh.14", "unc.42")
    yy = proportions;
    yy = rbind(rep(0, ncol(proportions)), proportions)
    yy = cbind(rep(1, nrow(yy)), yy)
    colnames(yy) = c("C0", colnames(proportions))
    rownames(yy) = c("background", rownames(proportions))
    proportions = yy
    
    save(expression, proportions, file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                                                version.Fraction.Matrix,"_", version.ExprsMatrix, ".Rdata"))
    
  }else{
    rownames(expression.sel) = c("ASE", "Serotonergic", "Dopaminergic", "Glutamatergic", "Ciliatedsensory",
                                  "GABAergic", "Mechanosensory", "unc.86", "Pharyngeal", "Cholinergic", "ceh.14", "unc.42")
    
  }
  #index.sel = c(13, 2, 1, 3, 6, 5, 7, 8, 9, 4, 10, 11)
  #proportions.sel = proportions[index.sel, ]
}

###############################
# select genes of intereset and match sample orders
###############################
load(file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed.Rdata"))
load(file = paste0(RdataDir, "Enrichscores_Matrix_13samples_selected_and_all_genes_", version.EnrichscoreMatrix, ".Rdata"))
enriched.list = colnames(enrich.matrix.sel)
#enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)

mm = match((enriched.list), colnames(expression))
expression.sel = expression[, mm]
#expression.sel = log2(expression.sel)

sels = match(rownames(expression.sel), rownames(proportions))
proportions.sel = proportions[sels, ]

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
  
  if(fitting.space == 'linear') {logaxis = 'xy'; yy = t(log2(t(yy)/yy[which(rownames(yy)=='background'), ]));
  }else{logaxis = ''}
  
  
  pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=TRUE, 
           color = c("lightgray", "blue"), legend = FALSE)
  
  pheatmap(yy, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, 
           cluster_cols=TRUE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))
  
  ## double check the expression matrix

  par(mfrow= c(1:2))
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
x=as.matrix(proportions.sel)
y = as.matrix(expression.sel)

x = x >0 

Example2test = c("lsy-6", "mir-791", "mir-793", "mir-790", "mir-1821", "mir-83", "mir-124")
jj2test = match(Example2test, colnames(y))
y = y[, jj2test[which(!is.na(jj2test)==TRUE)]]

res = matrix(NA, nrow = ncol(x), ncol = ncol(y)) 
colnames(res) = colnames(y)
rownames(res) = colnames(x)
#x.ms = apply(x, 2, sum)
#x = x[, which(x.ms>0)]
#x = x>0
alpha = 0.7 # 0.001 is the default value

intercept=FALSE
standardize=FALSE ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
standardize.response=FALSE
grouped = FALSE

for(n in 1:ncol(y))
{
  # n = 2
  cv.fit=cv.glmnet(x, y[,n], family='gaussian', alpha=alpha, nlambda=200, standardize=standardize, lower.limits = 0,
                   standardize.response=standardize.response, intercept=intercept, grouped = FALSE)
  fit=glmnet(x,y[,n], alpha=alpha, lambda=cv.fit$lambda,family='gaussian', lower.limits = 0,
             standardize=standardize, standardize.response=standardize.response, intercept=intercept)
  
  #par(mfrow= c(1,1))
  #plot(cv.fit)
  #plot(fit, label = TRUE)
  #plot(fit, xvar = "lambda", label = TRUE); abline(v=log(cv.fit$lambda.min))
  
  #myCoefs <- coef(fit, s=cv.fit$lambda.min);
  myCoefs = coef(cv.fit, s="lambda.min")
  res[,n] = as.numeric(myCoefs)[-1]
  
  cat("---------------\n", colnames(res)[n], ":\n", 
      paste0(rownames(res)[which(res[,n]>0)], collapse = "\n"), "\n",
      "---------------\n")
    
  #coef(cv.fit, s = "lambda.min")
  #myCoefs[which(myCoefs != 0 ) ]               #coefficients: intercept included
  #myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ] #feature names: intercept included
}

pdfname = paste0(resDir, "/res_deconv", 
                 "_fitting.", fitting.space, "_alpha.", alpha,  version.analysis, ".pdf")
pdf(pdfname, width=15, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))

pheatmap(log2(res[-1, ]+2^-10), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
         cluster_cols=FALSE, main = paste0("alpha = ", alpha),
         color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))

dev.off()



######################################
######################################
## Section: Test known packages for deconvolution
######################################
######################################
Test.some.known.packages.deconvolution.methods = FALSE
if(Test.some.known.packages.deconvolution.methods){

  ## "dtangle" is R package to deconvolving cell type proportions from DNA microarray data. 
  ## https://github.com/gjhunt/dtangle/blob/master/vign/basic-deconvolution.md
  library('dtangle')
  names(shen_orr_ex)
  #library('dtangle.data')
  #names(shen_orr_ex)
  truth = shen_orr_ex$annotation$mixture
  pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
  })
  names(pure_samples) = colnames(truth)
  pure_samples
  
  Y <- shen_orr_ex$data$log
  Y[1:4,1:4]
  
  marker_list = find_markers(Y,pure_samples,data_type="microarray-gene",marker_method='ratio')
  
  lapply(marker_list$L,head)
  
  q = .1
  quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
  K = length(pure_samples)
  n_choose = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
  n_choose
  
  marks = marker_list$L
  
  dc <- dtangle(Y,pure_samples,n_choose,data_type='microarray-gene',markers=marks)
  
  phats <- dc$estimates
  plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
  abline(coef=c(0,1))
  
  #### test CellMix and DSA
  require(CellMix)
  acr <- ExpressionMix("GSE20300", verbose = 2)
  #x <- ExpressionMix('GSE19830', verbose=TRUE)
  #annotation(x)
  acr
  #res <- gedBlood(acr, verbose = TRUE)
  s <- esApply(acr, 1L, sd, na.rm = TRUE)
  i <- order(s, decreasing = TRUE)[1:5000]
  rescs <- ged(acr[i, ], coef(acr), method = "csSAM", data = acr$Status, nperms = 200,
               verbose = TRUE)
  x <- rmix(3, 100, 20, markers = 5)
  basismap(x, Rowv = NA)
  p0 <- abs(coef(x) + rmatrix(coef(x), dist = rnorm, sd = 0.2))
  # rescale into proportion (columns sum up to one)
  p0 <- scoef(p0)
  
  profplot(x, p0)
  # fit DSection MCMC model
  ds <- ged(x, p0, "DSection", maxIter = 20, verbose = TRUE)
  
  # extract mixed samples
  mix <- mixedSamples(x)
  # load TIGER marker list
  ml <- MarkerList('TIGER')
  ml
  
  names(ml)
  basisnames(x)
  
  ml <- ml[c('brain', 'liver', 'lung')]
  summary(ml)
  
  mlx <- convertIDs(ml, mix, verbose=TRUE)
  summary(mlx)
  
  profplot(mlx[,1:10], mix)
  
  mlsc <- extractMarkers(mlx, expb(mix, 2), method='SCOREM', alpha=10^-12)
  summary(mlsc)
  
  
}





