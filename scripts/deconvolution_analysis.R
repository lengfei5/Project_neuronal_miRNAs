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

fitting.space = "log2" ## linear or log2 transformed for expression matrix
Use.coarse.neuronClass.FractionMatrix = FALSE

Data.complete = TRUE
Use.mergedFractionMatrix = TRUE
add.background.sample.in.fitting.linear.space = TRUE
Use.mergedExpressionMatrix = FALSE # group the genes if they show similar gene expression pattern

Manually.unifiy.sample.names.forMatrix = TRUE
Check.ProprotionMatrix.ExpressionMatrix = FALSE

version.ExprsMatrix = "miRNAs_neurons_v1_2018_03_07"
version.Fraction.Matrix = "_miRNAs_neurons_20180525"
version.EnrichscoreMatrix = "20180506"

version.analysis = "_20190107"

RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")
resDir = "../results/decomvolution_results/"
if(!dir.exists(resDir)) dir.create(resDir)

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
  if(!Use.coarse.neuronClass.FractionMatrix){
    load(file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix", 
                       version.Fraction.Matrix, ".Rdata"))
    #write.csv(proportions, file = "../results/tables_for_decomvolution/proportion_matrix_for_Chiara.csv")
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
  }else{
    
    load(file = paste0(RdataDir, "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix", 
                       version.Fraction.Matrix, ".Rdata"))
    source('miRNAseq_functions.R')
    library("openxlsx")
    coarseClass = read.xlsx("../data/proportion_matrix_for_Coarse_neuron_classes.xlsx", rowNames = TRUE,
                       sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, skipEmptyCols = TRUE, na.strings = "NA")
    mm = match(c("SENSORY", "MOTOR", "INTER"), rownames(coarseClass))
    coarseClass = as.matrix(coarseClass[mm, ])
    
    coarseClass[which(is.na(coarseClass))] = 0
    coarseClass[which(coarseClass== "1?")] = 1;
    
    coarseClass = data.frame(coarseClass, stringsAsFactors = TRUE)
    cc = t(coarseClass)
    cc = matrix(as.numeric(cc), nrow = nrow(cc), ncol = ncol(cc))
    colnames(cc) = rownames(coarseClass)
    rownames(cc) = colnames(coarseClass)
    
    TEST.newcc = FALSE
    if(TEST.newcc){
      n = 14;
      j = 2;
      kk = which(cc[, j]>0)
      sum(proportions[n,kk])
    }
    
    newcc = proportions %*% cc
    
    rownames(newcc) = c("Dopaminergic","Serotonergic","Glutamatergic", "Cholinergic", "GABAergic",  "Ciliatedsensory",
                              "Mechanosensory", "unc.86", "Pharyngeal", "ceh.14", "unc.42", "unc.3","ASE", "Pan.neurons")
    
    save(newcc, file = paste0(RdataDir,
                              "Tables_Coarse_neuronClasses_FractionMatrix_for_Sensory_Motor_Inter", 
                              version.Fraction.Matrix, ".Rdata"))
  }
  
}

######################################
######################################
## Section: enrichment score matrix and expression matrix
######################################
######################################
if(!Use.mergedExpressionMatrix){
  
  #load(file = paste0(RdataDir, 'piRANormalized_cpm.piRNA_batchCorrectedCombat_reAveraged_', version.ExprsMatrix, '.Rdata'))
  load(file = paste0(RdataDir, 'piRANormalized_cpm.piRNAnorm_batchCorrectedCombat_calibratedProtEff_', version.ExprsMatrix, '.Rdata'))
  
  if(fitting.space == "log2"){
    cpm.piRNA.bc = log2(cpm.piRNA.bc.prot)
    cpm.piRNA.bc.meanrep = average.biological.replicates(cpm.piRNA.bc)
    jj = grep('_untreated', colnames(cpm.piRNA.bc.meanrep))
    total = apply(cpm.piRNA.bc.meanrep[, jj], 1, median)
    xx = data.frame(total, cpm.piRNA.bc.meanrep[, -jj])
    ncs = sapply(colnames(xx)[-c(1:2)], function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
    ncs = sapply(ncs, function(x) gsub("*.neurons", "", x), USE.NAMES = FALSE)
    colnames(xx) = c('whole.body', 'background', ncs)
    expression = xx[, -c(1)] # ignore the gene expression in the whole body
    
  }else{
    source('miRNAseq_functions.R')
    cpm.piRNA.bc.meanrep = average.biological.replicates(cpm.piRNA.bc.prot)
    jj = grep('_untreated', colnames(cpm.piRNA.bc.meanrep))
    total = apply(cpm.piRNA.bc.meanrep[, jj], 1, median)
    xx = data.frame(total, cpm.piRNA.bc.meanrep[, -jj])
    ncs = sapply(colnames(xx)[-c(1:2)], function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
    ncs = sapply(ncs, function(x) gsub("*.neurons", "", x), USE.NAMES = FALSE)
    colnames(xx) = c('whole.body', 'background', ncs)
    
    expression = xx[, -c(1)]
    
    expr.vars = variance.biological.replicates(cpm.piRNA.bc.prot)
    xx = data.frame(expr.vars[, -1])
    ncs = sapply(colnames(xx)[-c(1)], function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
    ncs = sapply(ncs, function(x) gsub("*.neurons", "", x), USE.NAMES = FALSE)
    colnames(xx) = c('background', ncs)
    
    expr.vars = xx;
    
    plot(as.vector(as.matrix(expression)), as.vector(as.matrix(expr.vars)), log = "xy")
    points(seq(1, 10^6, by = 100), seq(1, 10^6, by = 100)^1.775*exp(-2.0), type = "l", col= 'red')
    
  }
  
  #write.csv(cpm.piRNA.bc, file = "/Volumes/groups/cochella/Chiara/table_normalized_piRNA_batchCorrected_replicates.csv")
  #write.csv(cpm.piRNA.bc.meanrep, file = "/Volumes/groups/cochella/Chiara/table_normalized_piRNA_batchCorrected_averageRepliciates.csv")
  ####################
  ## here we transform the gene expression by e'= (expression-background)/background 
  # then e' = 0 if e'<0 or e'<1; 
  # this transformation is to scale the range for each gene so that all data fall into the same range and e' should still follow the positive constrain 
  # meanwhile using ratio between expression and background to filter non-expressed ones  
  ###################
  
  Save.table.for.Chiara = FALSE
  if(Save.table.for.Chiara){
    write.csv(cpm.piRNA.bc.prot, 
              file = paste0("../results/final_tables_4Chiara/", 
                            "expression_matrix_allReplicates_piRNANorm_bacthCorrection_promoterCalibration", 
                            version.analysis, ".csv"), col.names = TRUE, row.names = TRUE)
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
 
  #old.sampleNames = rownames()
  newSampleNames.prop = c("Dopaminergic","Serotonergic","Glutamatergic", "Cholinergic", "GABAergic",  "Ciliatedsensory",
                           "Mechanosensory", "unc.86", "Pharyngeal", "ceh.14", "unc.42", "unc.3","ASE", "Pan.neurons")
  print(cbind(rownames(proportions), newSampleNames.prop))
  rownames(proportions) = newSampleNames.prop;
  
  expression = t(expression)
  newSampleNames.expr = c( "background", "ASE", "Serotonergic", "Dopaminergic", "Glutamatergic", "Ciliatedsensory",
                           "GABAergic", "Mechanosensory", "unc.86", "Pharyngeal", "Cholinergic", "ceh.14", "unc.42", 
                           "unc.3", "Pan.neurons")
  print(cbind(rownames(expression), newSampleNames.expr))
  rownames(expression) = newSampleNames.expr
  
  ###############################
  # if fitting in linear scale, then the background sample is used and the background is added in proportion matrix as well
  # if fitting in log2 scale, then the log2(sample/background) is used; in this case, the coefficients must be also positive
  ###############################
  if(fitting.space == "linear"){
    
    expr.vars = t(expr.vars)
    print(cbind(rownames(expr.vars), newSampleNames.expr))
    rownames(expr.vars) = newSampleNames.expr
    
    if(add.background.sample.in.fitting.linear.space){
      yy = proportions;
      yy = rbind(rep(0, ncol(proportions)), proportions)
      #yy = cbind(rep(1, nrow(yy)), yy)
      #colnames(yy) = c("C0", colnames(proportions))
      rownames(yy) = c("background", rownames(proportions))
      ss = apply(yy, 2, sum)
      yy = yy[, which(ss>0)]
      
      proportions = yy
    }
    
    save(expression, proportions, expr.vars, file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                                                version.Fraction.Matrix,"_", version.ExprsMatrix, "fitting_scale", 
                                                fitting.space,".Rdata"))
    
  }else{
    yy = expression;
    # log2 (expresion/background)
    #for(n in 1:nrow(yy)) yy[n, ] = log2(as.numeric(yy[n, ])/as.numeric(expression[1,]))
    for(n in 1:nrow(yy)) yy[n, ] = as.numeric(yy[n, ] - as.numeric(expression[1,]))
    
    yy = yy[-1, ]
    expression = yy;
    
    save(expression, proportions, file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                                                version.Fraction.Matrix,"_", version.ExprsMatrix, "fitting_scale", 
                                                fitting.space,".Rdata"))
  }
  
}

########################################################
########################################################
# Section :
# select genes of intereset and match sample orders
########################################################
########################################################
load(file = paste0(RdataDir, "Expression_Fraction_Matrix_withBackground_cleaed", 
                     version.Fraction.Matrix,"_", version.ExprsMatrix, "fitting_scale", 
                     fitting.space,".Rdata"))
exprThreshold.pan.vs.bg.log2 = 2.0

# use coarse neuron group or not
if(Use.coarse.neuronClass.FractionMatrix){
  load(file = paste0(RdataDir, "Tables_Coarse_neuronClasses_FractionMatrix_for_Sensory_Motor_Inter", 
                version.Fraction.Matrix, ".Rdata"))
  proportions = newcc;
}

### select genes of interest
load(file = paste0(RdataDir, "Enrichscores_Matrix_13samples_selected_and_all_genes_", version.EnrichscoreMatrix, ".Rdata"))
enriched.list = colnames(enrich.matrix.sel)
#enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)

mm = match((enriched.list), colnames(expression))
expression.sel = expression[, mm]

if(fitting.space == "log2"){
  sel.pan.neurons = which(expression.sel[which(rownames(expression.sel)=="Pan.neurons"), ] > exprThreshold.pan.vs.bg.log2)
  expression.sel = expression.sel[, sel.pan.neurons]
}else{
  pans = expression.sel[which(rownames(expression.sel)=="Pan.neurons"), ];
  bgs = expression.sel[which(rownames(expression.sel)=="background"), ];
  rrs = log2(pans/bgs)
  rrs[order(rrs)]
  
  sel.pan.neurons = which(log2(pans/bgs) > exprThreshold.pan.vs.bg.log2)
  
  expression.sel = expression.sel[, sel.pan.neurons]
}

sels = match(rownames(expression.sel), rownames(proportions))
proportions.sel = proportions[sels, ]

if(Data.complete){
  jj = match(rownames(proportions), rownames(proportions.sel))
  proportions.sel = proportions.sel[jj, ]
  expression.sel = expression.sel[jj, ]
}

####################
## double check the proprotion matrix and expression matrix
####################
if(Check.ProprotionMatrix.ExpressionMatrix){
  source("miRNAseq_functions.R")
  
  pdfname = paste0(resDir, "/Heatmap_Proportiona_Expression_Matrix", version.analysis, ".pdf")
  pdf(pdfname, width=15, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  Plot.ProprotionMatrix.ExpressionMatrix(proportions.sel, expression.sel, fitting.space = "log2")
    
  dev.off()
  
  Save.table.for.Chiara = FALSE
  if(Save.table.for.Chiara){
    write.csv(proportions.sel, 
              file = paste0("../results/final_tables_4Chiara/", 
                            "proportionas_matrix_used4deconvolution", version.analysis, ".csv"), col.names = TRUE, row.names = TRUE)
    write.csv(expression.sel, 
              file = paste0("../results/final_tables_4Chiara/", 
                            "expression_matrix_used4deconvolution", version.analysis, ".csv"), col.names = TRUE, row.names = TRUE)
    
    kk = match(colnames(expression), colnames(expression.sel))
    write.csv(expression[, which(is.na(kk))], 
              file = paste0("../results/final_tables_4Chiara/", 
                            "expression_matrix_NOTused4deconvolution_notEnriched_and_lowlyExpressed", 
                            version.analysis, ".csv"), col.names = TRUE, row.names = TRUE)
    
  }
 
}

######################################
######################################
# Section: process the table and run the glmnet
# the glmnet will be run for each gene, because the group lasso is not desirable
######################################
######################################
x = as.matrix(proportions.sel)
y = as.matrix(expression.sel)

if(fitting.space == "log2"){
  ## force value to be zero if they are lower than the background
  y[y<0] = 0
}

if(!Use.coarse.neuronClass.FractionMatrix){
  
  x = x > 0
  ss = apply(x, 2, sum)
  x = x[, which(ss>0)]
  Example2test = c("lsy-6", "mir-791", "mir-793",  "mir-792","mir-1821", "mir-83", "mir-124")
  #Example2test = c(Example2test, setdiff(colnames(y), Example2test))
  jj2test = match(Example2test, colnames(y))
  y = y[, jj2test[which(!is.na(jj2test)==TRUE)]]
  
  if(fitting.space == "linear"){
    weights = expr.vars[, match(colnames(y), colnames(expr.vars))]
    weights = weights[match(rownames(y), rownames(weights)), ]
    weights = (1/weights)^0.5
  }
  
  res = matrix(NA, nrow = ncol(x), ncol = ncol(y)) 
  colnames(res) = colnames(y)
  rownames(res) = colnames(x)
  
  ##########################################
  # glmnet with global alpha parameter or gene-specific alpha parameters
  ##########################################
  library("pheatmap")
  library("RColorBrewer")
  TEST.glmnet.gene.specific.alpha = FALSE
  save.deconvolution.results.for.downstream.analysis = TRUE
  
  Regroup.highly.correlated.neurons = TRUE
  Test.groupLasso = FALSE;
  
  if(Regroup.highly.correlated.neurons){
    source("select_tuningParams_elasticNet.R")
    x = Regroup.Highly.Correlated.Neurons(x, cor.cutoff = 0.8, plot.grouping.result = FALSE)
  }
  #Methods2test = c("cv.lambda.1se", "cv.lambda.min", "bic", "aic", "aicc")
  #Methods2test = c("cv.lambda.1se", "bic")
  Methods2test = c("cv.lambda.1se")
  #alphas = c(seq(0.1, 1, by= 0.1))
  alphas = c(0.005, seq(0.01, 0.1, by= 0.01), seq(0.2, 1, by=0.1))
  #alphas = c(0.1)
  lambda = 10^seq(-3, 3, length.out = 500)
  nlambda = 500;
  
  # make a folder for the result
  if(TEST.glmnet.gene.specific.alpha) {
    alpha.hyperparam = "gene.specific.alpha"
  }else{
    alpha.hyperparam = "global.alpha"
  }
  
  testDir = paste0(resDir, "deconv_results")
  
  if(!dir.exists(testDir)) system(paste0('mkdir -p ', testDir))
  
  
  source("select_tuningParams_elasticNet.R")
  for(method in Methods2test)
  {
    cat("-- model selection method -- ", method, "\n")
    pdfname = paste0(testDir, "/deconv_res", 
                     "_fitting.", fitting.space, 
                     "_glmnet_global_alpha_method_select_tuning_parameters_", method, "_", alpha.hyperparam, version.analysis, 
                     ".pdf")
    
    pdf(pdfname, width=22, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    par(mfrow=c(1, 1))
    
    if(!Test.groupLasso){
      keep = run.glmnet.select.tuning.parameters(x, y, alphas = alphas, method = method, lambda = lambda, intercept = TRUE, standardize = TRUE, nfold = 7, 
                                                 Gene.Specific.Alpha = TEST.glmnet.gene.specific.alpha);
    }else{
      #source("select_tuningParams_elasticNet.R")
      keep = run.gglasso.select.tuning.parameters(x, y, cor.cutoff=seq(1, 0.5, by= -0.1), method = method, lambda = lambda, intercept = TRUE, nfold = 7)
      
    }
    
    dev.off()
    
    save(x, y, alphas, keep, file = paste0(RdataDir, "deconvolution_results_glmnet_log2scale_method_", method, "_", alpha.hyperparam,
                                           version.analysis, ".Rdata"))
    
  }
  
  
  if(save.deconvolution.results.for.downstream.analysis){
    
    method = "cv.lambda.1se"
    alpha.hyperparam = "global.alpha"
    load(file = paste0(RdataDir, "deconvolution_results_glmnet_log2scale_method_", method, "_", alpha.hyperparam,
                       version.analysis, ".Rdata"))
    
    tabDir = paste0(resDir, "deconv_results/tables_logscale/")
    if(!dir.exists(tabDir)) system(paste0('mkdir -p ', tabDir))
    
    pdfname = paste0(tabDir, "deconvolution_results_glmnet_log2scale_method_", 
                     method, "_", alpha.hyperparam,
                     version.analysis, ".pdf")
    pdf(pdfname, width=20, height = 12)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    par(mfrow=c(1, 1))
    
    for(n in 1:length(alphas)){
      # n = 9;
      cat("alpha -- ", alphas[n], "\n")
      res = keep[[n]];
      
      source("select_tuningParams_elasticNet.R")
      res = clustering.gene.neuronClass(res);
      
      write.csv(res, file = paste0(tabDir, "deconvolution_results_glmnet_log2scale_method_", method, "_", alpha.hyperparam, 
                                   version.analysis, "_alpha_", alphas[n], ".csv"), row.names = TRUE)
      
      cols = colorRampPalette((brewer.pal(n = 7, name="Reds")))(100)
      pheatmap(res, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE,
               cluster_cols=FALSE, main = paste0("alpha = ", alphas[n], " -- ", method), na_col = "white",
               color = cols)
      
    }
    
    dev.off()
    
  }
  
}else{
  # write.csv(x, file = "/Volumes/groups/cochella/Chiara/table_cellNbs_in_Sensory_Motor_Inter_for_14_Samples.csv")
  coarseGroup_Dir = paste0(resDir, "deconv_results_log2_coarse")
  if(!dir.exists(coarseGroup_Dir)) system(paste0('mkdir -p ', coarseGroup_Dir))
  
  pdfname = paste0(coarseGroup_Dir, "/deconv_coaseGroup",
                   "_fitting.", fitting.space, 
                   "_Lasso",  version.analysis, ".pdf")
  
  pdf(pdfname, width=20, height = 12)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  par(mfrow=c(1, 1))
  
  source("select_tuningParams_elasticNet.R")
  keep.coarse = run.glmnet.for.coarse.groups(x, y)
  
  
  res[which(is.na(res))] = 0 
  cols = c("white", colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
  pheatmap(res, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=TRUE, main = paste0(" coarse group of neurons"), na_col = "white",
           color = cols)
  
  dev.off()
  
  save(keep.coarse, file = paste0(RdataDir, "preliminary_results_for_Coarse_neuronGroups", version.analysis, ".Rdata"))
  
  res = keep.coarse;
 
  
  write.csv(res, file = paste0(coarseGroup_Dir, "/deconvolution_results_glmnet_Lasso_method_fitting_", fitting.space, 
                               "_", version.analysis, ".csv"), row.names = TRUE)
  
}



