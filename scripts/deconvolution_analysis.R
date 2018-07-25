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

resDir = "../results/decomvolution_results"
if(!dir.exists(resDir)) dir.create(resDir)

fitting.space = "linear" ## linear or log2 transformed for expression matrix
Use.mergedFractionMatrix = TRUE
Use.mergedExpressionMatrix = FALSE
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
  load(file = paste0(RdataDir, "Enrichscores_Matrix_13samples_selected_and_all_genes_", version.EnrichscoreMatrix, ".Rdata"))
  enriched.list = colnames(enrich.matrix.sel)
  #enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)
  
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
  expression = xx[, -c(1)]
  
  for(n in 1:ncol(expression)) 
  {
    if(fitting.space == "linear"){
      expression[,n] = (expression[,n]-xx$background)/xx$background
      ## use the e'> 1 as a threshold 
      expression[which(expression[,n]<1),n] = 0
      #expression[which(expression<1)] = 0
    }else{
      expression[,n] = log2(expression[,n]/xx$background)
    }
  }
  
}else{
  source('miRNAseq_functions.R')
  newExprM = expressionMatrix.grouping(xx) 
}

####################
## double check the proprotion matrix and expression matrix
## match the sample order in the proprotion matrix and expression matrix 
## now manually (to change)
####################

####################
## select the miRNAs to analyze and match the samples for fraction matrix and expression matrix 
####################
mm = match((enriched.list), rownames(expression))
expression.sel = t(expression[mm, ])
#expression.sel = log2(expression.sel)

index.sel = c(13, 2, 1, 3, 6, 5, 7, 8, 9, 4, 10, 11)
proportions.sel = proportions[index.sel, ]

if(Check.ProprotionMatrix.ExpressionMatrix){
  xx = proportions;
  xx[which(xx>0)] = 1
  yy = expression.sel[c(3, 2, 4, 10, 6, 5, 7, 8, 9, 11, 12, 1), ];
  #yy = proportions[c(index.sel, 12, 14), ]
  
  library("pheatmap")
  library("RColorBrewer")
  
  pdfname = paste0(resDir, "/heatmap_for_ProportionaMatrix_14_samples_96_neuronClasses_ExpressionMatrix", ".pdf")
  pdf(pdfname, width=15, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  # par(mfcol=c(1, 1))
  
  pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=TRUE, 
           color = c("lightgray", "blue"), legend = FALSE)
  
  pheatmap(yy, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, 
           cluster_cols=TRUE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))
  
  ## double check the expression matrix
  if(fitting.space == 'linear') {logaxis = 'xy';
  }else{logaxis = ''}
  
  par(mfrow= c(1:2))
  plot(t(expression.sel[match(c("Dopaminergic", "Ciliated.sensory"), rownames(expression.sel)), ]), log=logaxis)
  abline(0, 1, lwd=2.0, col='red')
  plot(t(expression.sel[match(c("mechanosensory",  "unc.86.expressing"), rownames(expression.sel)), ]), log=logaxis)
  abline(0, 1, lwd=2.0, col='red')
  
  dev.off()
  
}


######################################
######################################
## Section: process the table and run the glmnet
######################################
######################################
require(glmnet)
x=as.matrix(proportions.sel)
y = as.matrix(expression.sel)

x.ms = apply(x, 2, sum)
x = x[, which(x.ms>0)]

#x = x>0
#y=cbind(y1,y2)
#rownames(y) = res.sel$gene

intercept=0
standardize=FALSE ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
standardize.response=FALSE
alpha = 0.5
grouped = FALSE

cv.fit=cv.glmnet(x, y, family='mgaussian', grouped=grouped, alpha=alpha, nlambda=500, standardize=standardize, 
                 standardize.response=standardize.response, intercept=intercept)
par(mfrow= c(1,1))
plot(cv.fit)

#optimal = which(cv.fit$lambda==cv.fit$lambda.min)
#optimal = which(cv.fit$lambda==cv.fit$lambda.1se)
fit=glmnet(x,y, alpha=alpha, lambda=cv.fit$lambda,family='mgaussian', 
            standardize=standardize, standardize.response=standardize.response, intercept=intercept)

#optimal = which(fit$df<=70)
#optimal = max(optimal)
optimal = which(cv.fit$lambda==cv.fit$lambda.min)

## Adapted from @Mehrad Mahmoudian:
myCoefs <- coef(fit, s=cv.fit$lambda.min);
#myCoefs[which(myCoefs != 0 ) ]               #coefficients: intercept included
## [1]  1.4945869 -0.6907010 -0.7578129 -1.1451275 -0.7494350 -0.3418030 -0.8012926 -0.6597648 -0.5555719
## [10] -1.1269725 -0.4375461
#myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ] #feature names: intercept included
## [1] "(Intercept)" "feature1"    "feature2"    "feature3"    "feature4"    "feature5"    "feature6"   
## [8] "feature7"    "feature8"    "feature9"    "feature10"  

## Asseble into a data.frame
#myResults <- data.frame(
#  features = myCoefs@Dimnames[[1]], #intercept included
#  coefs    = myCoefs              #intercept included
#)

#myResults
##       features      coefs

## extract the fitting result
#colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
res = matrix(NA, nrow = ncol(x), ncol = ncol(y))
colnames(res) = colnames(y)
rownames(res) = colnames(x)
for(n in 1:ncol(res))
{
  res[,n] = fit$beta[[n]][,optimal]
}

res = data.frame(res)

cbind(rownames(res)[order(-res$lsy.6)],  res$lsy.6[order(-res$lsy.6)])



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





