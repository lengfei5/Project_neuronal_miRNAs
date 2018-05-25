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
RdataDir = paste0("../results/miRNAs_neurons_v1_2018_03_07/Rdata/")

version.ExprsMatrix = "miRNAs_neurons_v1_2018_03_07"

resDir = "../results/decomvolution_results"
if(!dir.exists(resDir)) dir.create(resDir)

fitting.space = "linear" ## linear or log2 transformed for expression matrix


######################################
######################################
## Section: import the fraction matrix and expression matrix
######################################
######################################


######################################
######################################
## Section: process the table and run the glmnet
######################################
######################################



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





