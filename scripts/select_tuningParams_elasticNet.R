##########################################################################
##########################################################################
# Project:
# Script purpose: to select the optimal model using BIC, AIC, AICc and HQC for elastic-net
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 19 11:09:08 2018
##########################################################################
##########################################################################
###############################
# function to estimate degree of freedom for elastic-net 
# this function is modified from 
# https://github.com/lcallot/lassovar/blob/master/R/lassovar-ada.R
# and the original paper comes from is following:
# DEGREES OF FREEDOM IN LASSO PROBLEMS By Ryan J. Tibshirani and Jonathan Taylor1
###############################
estimate.elasticNet.df <-function(xx, alpha = 0,  fit)
{
  xx <- (cbind(1, xx))
  dfs <- NULL
  
  for(l in fit$lambda){
    # l = fit$lambda[1]
    ccf = coef(fit, s = l)
    kk = which(ccf != 0)
    if(length(kk)==0) {
      dfs = c(dfs, 0)
    }else{
      if(length(kk)==1) {
        dfs = c(dfs, 1)
      }else{
        dfs <-c(dfs, sum(diag(xx[,kk] %*% solve(t(xx[ ,kk]) %*% xx[,kk] + l*(1.0-alpha)*diag(ncol(xx[,kk])))%*% t(xx[,kk]))))
      }
    }
  }
  
  #plot(fit$lambda, dfs, col='blue')
  #points(fit$lambda, fit$df, col='red')
  
  #abline(0, 1, lwd=2.0, col = 'red')
  return(dfs)
}

calcluate.loglikelihood.glmnet = function(xx, yy, fit, use.deviance.glmnet = TRUE)
{
  if(use.deviance.glmnet){
    loglike = -0.5*deviance.glmnet(fit)
  }else{
    model = fit;
    coef = coef(model)
    lambda = model$lambda
    #df = model$df
    
    yhat=cbind(1,xx)%*%coef
    
    residuals = (yy- yhat)
    mse = colMeans(residuals^2)
    #sse = colSums(residuals^2)
    loglike = log(mse)*length(yy)*(-0.5)
  }
  
  return(loglike)
}

###############################
# funtion to select tuning parameter for glmnet using methods, 
# cv.lambda.1se, bic, aic, ebic, hbic ...
###############################
select.tuning.parameters.for.glmnet = function(xx, yy, alpha = 0.1, fit, cv.fit, method = c("cv.lambda.1se"))
{
  # xx = x; yy = y[,n];
  # coef(fit, s=cv.fit$lambda.1se);
  if(method == "cv.lambda.1se") {
    myCoefs <- coef(fit, s=cv.fit$lambda.1se);
  }else{
    if(method == "cv.lambda.min") {
      myCoefs <- coef(fit, s=cv.fit$lambda.min);
      
    }else{
      #n.var = (fit$df + 1);
      n.data = length(yy);
      n.var = estimate.elasticNet.df(xx, alpha = alpha,  fit)
      # plot((fit$df+1), n.var); abline(0, 1, lwd=2.0, col='red')
      #loglike.glment = calcluate.loglikelihood.glmnet(xx, yy, fit, use.deviance.glmnet = TRUE)
      loglike = calcluate.loglikelihood.glmnet(xx, yy, fit, use.deviance.glmnet = FALSE)
      
      # plot(loglike.glment, loglike.mse);abline(0, 1, lwd=2.0, col='red')
      
      if(method == "bic"){
        crits = -2*loglike + n.var *log(n.data)
        optim.crits = which(crits==min(crits))
      }
      if(method == "aic"){
        crits = -2*loglike + 2*n.var;
        optim.crits = which(crits==min(crits))
      }
      if(method == "aicc"){
        crits = -2*loglike + 2*n.var + (2*nvar*(nvar+1))/(n.data-nvar-1);
        optim.crits = which(crits==min(crits))
      }
      
      #plot(fit, label = TRUE); abline(v= c(cv.fit$lambda.min, cv.fit$lambda.1se)); abline(v=c(fit$lambda[optim.crits]), col='red')
      myCoefs <- coef(fit, s=fit$lambda[optim.crits]);
      
    }
  }
  
  #beta[n, (ncol(x)+1)] = mse
  #beta[n, (ncol(x)+2)] =  cv.fit$cvm[n] # cv
  #beta[n, (ncol(x)+3)] = n.data*log(mse) + 2*nvar # aic
  #beta[n, (ncol(x)+4)] = n.data*log(mse) + 2*nvar + (2*nvar*(nvar+1))/(n.data-nvar-1) # aicc
  #beta[n, (ncol(x)+5)] = n.data*log(mse) + 2*nvar*log(log(n.data)) # hqc
  #beta[n, (ncol(x)+6)] = n.data*log(mse) + nvar*log(n.data) # bic
  #gamma = 0.5; beta[n, (ncol(x)+7)] = n.data*log(mse) + nvar*(log(n.data) + 2*gamma*log(p_full)) # ebic
  #gamma = 0.5; beta[n, (ncol(x)+8)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic
  #gamma = 1; beta[n, (ncol(x)+9)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic2
  #gamma = 2; beta[n, (ncol(x)+10)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic3
  
  return(myCoefs)
  
}

########################################################
########################################################
# Section: test glmnet with global alpha and gene-specific alpha
########################################################
########################################################
run.glmnet.select.tuning.parameters = function(x, y, alphas = seq(0.1, 0.5, by=0.1), lambda = NULL, nfold = 10, 
                                               nlambda = 100, intercept = TRUE, standardize = FALSE,
                                               Gene.Specific.Alpha = FALSE,  
                                               method = "cv.lambda.1se",   
                                               plot.cluster.col = FALSE,
                                               omit.BIC = TRUE,  plot.crit.bic = TRUE)
{
  require(glmnet)
  library("pheatmap")
  library("RColorBrewer")
  
  # basic arguments for cv.glmnet and glmnet
  #intercept=TRUE
  #standardize=TRUE ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
  standardize.response=FALSE
  grouped = FALSE
  cols = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100)
  
  if(!Gene.Specific.Alpha){
    ###############################
    # global alpha parameter is tested in glmnet
    # so start the loop with alphas 
    # but it does not work well
    ###############################
    results = list();
    #ii.alpha = 0;
    for(alpha in alphas)
    {
      cat("global alpha --", alpha, "----------\n")
      # alpha = 0.7; lambda = NULL; nlambda = 200; method = "bic"
      
      res = matrix(NA, nrow = ncol(x), ncol = ncol(y)) 
      colnames(res) = colnames(y)
      rownames(res) = colnames(x)
      set.seed(1)
      for(n in 1:ncol(y))
      {
        # n = 1
        if(is.null(lambda)){
          cv.fit=cv.glmnet(x, y[,n], family='gaussian', alpha=alpha, nlambda=nlambda, standardize=standardize, lower.limits = 0,
                           standardize.response=standardize.response, intercept=intercept, grouped = FALSE) 
          
        }else{
          cv.fit=cv.glmnet(x, y[,n], family='gaussian', alpha=alpha, lambda = lambda, standardize=standardize, lower.limits = 0,
                           standardize.response=standardize.response, intercept=intercept, grouped = FALSE)
        }
        
        #par(mfrow= c(1,1))
        # plot(cv.fit, main = colnames(y)[n])
        
        fit=glmnet(x,y[,n], alpha=alpha, lambda=cv.fit$lambda,family='gaussian', lower.limits = 0,
                   standardize=standardize, standardize.response=standardize.response, intercept=intercept)
        
        myCoefs = select.tuning.parameters.for.glmnet(x, y[,n], alpha, fit, cv.fit, method = method)
        
        res[, n] = as.numeric(myCoefs)[-1]
        
        #cat("---------------\n", colnames(res)[n], ":\n", 
        #    paste0(rownames(res)[which(res[,n]>0)], collapse = "\n"), "\n",
        #    "---------------\n")
      }
      
      if( ! plot.cluster.col){
        res[which(res==0)] = NA;
        pheatmap((res), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
                 cluster_cols=plot.cluster.col, main = paste0("alpha = ", alpha, " - method : ", method), na_col = "gray30",
                 color = cols)
      }else{
        pheatmap((res), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
                 cluster_cols=plot.cluster.col, main = paste0("alpha = ", alpha, " - method : ", method), na_col = "gray30",
                 color = cols)
        res[which(res==0)] = NA;
      }
      
      results[[which(alphas==alpha)]] = res; 
    }
    
    return(results)
    
  }else{
    ###############################
    # select gene-specific alpha parameter
    # so the selection will be done gene by gene
    # the loop starts with gene
    # this can be done by home-made, or with examples 
    # https://www.r-bloggers.com/variable-selection-with-elastic-net/
    # or https://stats.stackexchange.com/questions/268885/tune-alpha-and-lambda-parameters-of-elastic-nets-in-an-optimal-way (CARET)
    # https://rpubs.com/ledongnhatnam/250381
    # https://daviddalpiaz.github.io/r4sl/elastic-net.html
    ###############################
    my.own.version = FALSE
    use.foreach = TRUE
    use.CARET = FALSE
    
    for(n in 1:ncol(y) ){
      # n = 1
      
      cat("start the grid search for tuning parameter alpha -----\n")
      cat("----n -- ", n, "--", colnames(y)[n],  " -- ")
      
      require(glmnet)
      library(doParallel)
      
      if(use.foreach){
        library(foreach)
        library(pROC)
        
        #pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
        #lapply(pkgs, require, character.only = T)
        registerDoParallel(cores = 6)
        
        set.seed(2017)
        #cv.fit=cv.glmnet(x, y[,n], family='gaussian', alpha=alpha, nlambda=nlambda, standardize=standardize, lower.limits = 0,
        #                 standardize.response=standardize.response, intercept=intercept, grouped = FALSE) 
        search <- foreach(i = alphas, .combine = rbind) %dopar% {
          cv <- cv.glmnet(x, y[,n], family = "gaussian", alpha = i, nfold = 10, nlambda=nlambda, standardize=standardize, lower.limits = 0,
                          standardize.response=standardize.response, intercept=intercept, grouped = FALSE, paralle = TRUE)
          data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
        }
        
        cv3 <- search[search$cvm == min(search$cvm), ]
        cat("selected alpha = ", cv3$alpha, " -- lambda = ", cv3$lambda.1se, "\n")
        
        fit=glmnet(x,y[,n], alpha=cv3$alpha, lambda=cv3$lambda.1se, family='gaussian', lower.limits = 0,
                   standardize=standardize, standardize.response=standardize.response, intercept=intercept)
        #fit <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
        #coef(fit)
        myCoefs <- coef(fit)
      }
      
      if(my.own.version){
        alphas = seq(0.1, 1.0, by = 0.1)
        lambda = 10^(seq(4, -2, length.out = 500))
        # nfolds = 5;
        #foldid=sample(1:nfolds,size=length(y[,n]),replace=TRUE)
        #lambdas.save = NULL;
        lambdas.min = c()
        lambdas.1se = c()
        cvms = c()
        
        for(m in 1:length(alphas))
        {
          #registerDoMC(cores=2)
          registerDoParallel(detectCores())
          set.seed(1010);
          
          alpha = alphas[m]
          #cvm = rep(0, length(lambda))
          
          cv.fit=cv.glmnet(x, y[,n], family='gaussian', alpha=alpha, lambda = lambda, nfolds = nfolds, 
                           standardize=standardize, lower.limits = 0, nlambda = length(lambda),
                           standardize.response=standardize.response, intercept=intercept, grouped = FALSE, parallel = TRUE)
          
          #for(k in 1:nb.cv){
          #  cat(length(cv.fit$cvm), "---\n")
          #  mm = match(cv.fit$lambda, lambda)
          #  cvm = cvm + cv.fit$cvm
          #}
          #cvm = cvm/nb.cv
          par(mfrow= c(1,1))
          plot(cv.fit, main = paste0(colnames(y)[n], "-- alpha : ", alpha))
          
          cvms = c(cvms, cv.fit$cvm[which(cv.fit$lambda==cv.fit$lambda.min)])
          #cvms = c(cvms, cv.fit$cvup[which(cv.fit$lambda==cv.fit$lambda.1se)])
          #lambdas[[m]] = cv.fit$lambda;
          lambdas.min = c(lambdas.min, cv.fit$lambda.min)
          lambdas.1se = c(lambdas.1se, cv.fit$lambda.1se)
          #plot(fit, label = TRUE)
          #plot(fit, xvar = "lambda", label = TRUE); abline(v=log(cv.fit$lambda.min))
        }
        
        index.alpha = which(cvms==min(cvms))
        fit=glmnet(x,y[,n], alpha=alphas[index.alpha], lambda=lambda, family='gaussian', lower.limits = 0,
                   standardize=standardize, standardize.response=standardize.response, intercept=intercept)
        #myCoefs <- coef(fit, s=lambdas.min[index.alpha]);
        myCoefs <- coef(fit, s=lambdas.1se[index.alpha]);
        
      }
      
      res[,n] = as.numeric(myCoefs)[-1]
      
      cat("---------------\n", colnames(res)[n], ":\n", 
          paste0(rownames(res)[which(res[,n]>0)], collapse = "\n"), "\n",
          "---------------\n")
      
    }
    
    res[which(res==0)] = NA;
    pheatmap((res), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
             cluster_cols=FALSE, main = paste0("gene-specific alpha - method : ", method), na_col = "gray30",
             color = cols)
    
    return(res)
  }
  
}

########################################################
########################################################
# Section : test gcdnet
# 
########################################################
########################################################
select.tuning.parameters.for.gcdnet = function(x, y,lambda = NULL, lambda2 = 0.1, method = "cv.lambda.1se", omit.BIC = FALSE,  
                                               nfold = 10, nlambda = 100, standardize = FALSE,  
                                               mse.overfitting.threshold = 10^-5, plot.crit.bic = TRUE, main = NULL)
{
  # y = y[,n] ; standardize = FALSE; nfold =10; nlambda =100;mse.overfitting.threshold = 10^-5
  if(is.null(lambda)){
    cv.fit=cv.gcdnet(x, y, nlambda=nlambda, lambda2 = lambda2, method = "ls", nfolds = nfold, standardize = standardize)
  }else{
    cv.fit=cv.gcdnet(x, y, lambda = lambda, lambda2 = lambda2, method = "ls", nfolds = nfold, standardize = standardize)
  }
  
  fit=gcdnet(x, y, lambda=cv.fit$lambda, method = "ls", lambda2 = lambda2, standardize=standardize)
  
  n.data = length(y)
  #n=length(y)
  list.lambda = fit$lambda;
  crits = c("mse", "cv", "AIC", "AICc", "HQC","BIC", "EBIC", "HBIC",  "HBIC2", "HBIC3")
  if(omit.BIC){
    
    if(method == "cv.lambda.min") {
      optimal = which(fit$lambda==cv.fit$lambda.min);
      coef = as.numeric(coef(fit, s = cv.fit$lambda.min))[-1]
    }
    if(method == "cv.lambda.1se") {
      optimal = which(fit$lambda==cv.fit$lambda.1se)
      coef = as.numeric(coef(fit, s = cv.fit$lambda.1se))[-1]
    }
    coef[which(coef==0)] = NA;
    
    result = list(coefficients = coef,
                  optimal = optimal)
    
  }else{
    beta = matrix(NA, nrow = length(list.lambda), ncol = (ncol(x)+length(crits)))
    #rownames(beta) = paste0(list.lambda);
    colnames(beta) = c(colnames(x), crits)
    
    for(n in 1:length(list.lambda))
    {
      # n = 40;
      #df = fit$df[n]
      cfs = as.numeric(coef(fit, s=list.lambda[n]))
      cfs = cfs[-1] # move the intercept
      feature.sels = which(cfs>0)
      names(feature.sels) = colnames(x)[feature.sels]
      
      if(length(feature.sels)==0) {
        mse = sum((y-mean(y))^2);
      }else{
        x.sel = x[, feature.sels]
        
        if(length(feature.sels)>=2){
          x.sel.y = data.frame(y, x.sel)
          model =  eval(parse(text = paste0("y ~ ", paste0(colnames(x.sel), collapse = " + "))))
          ff = lm(model, data = x.sel.y)
        }else{
          ff = lm(y ~ x.sel)
        }
        mse = sum(ff$residuals^2)/n.data
        beta[n, feature.sels] = cfs[feature.sels];
        #beta[n, feature.sels] = as.numeric(coefficients(ff))[-1]
      }
      
      nvar = fit$df[n] + 1
      p_full = ncol(x)
      
      beta[n, (ncol(x)+1)] = mse
      beta[n, (ncol(x)+2)] =  cv.fit$cvm[n] # cv
      beta[n, (ncol(x)+3)] = n.data*log(mse) + 2*nvar # aic
      beta[n, (ncol(x)+4)] = n.data*log(mse) + 2*nvar + (2*nvar*(nvar+1))/(n.data-nvar-1) # aicc
      beta[n, (ncol(x)+5)] = n.data*log(mse) + 2*nvar*log(log(n.data)) # hqc
      beta[n, (ncol(x)+6)] = n.data*log(mse) + nvar*log(n.data) # bic
      gamma = 0.5; beta[n, (ncol(x)+7)] = n.data*log(mse) + nvar*(log(n.data) + 2*gamma*log(p_full)) # ebic
      gamma = 0.5; beta[n, (ncol(x)+8)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic
      gamma = 1; beta[n, (ncol(x)+9)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic2
      gamma = 2; beta[n, (ncol(x)+10)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic3
    }
    
    beta = data.frame(beta)
    #plot(fit, label = TRUE)
    #plot(fit, xvar = "lambda", label = TRUE); abline(v=log(cv.fit$lambda.min))
    kk = which(beta$mse>mse.overfitting.threshold)
    if(length(kk)>0){
      beta.fitered = beta[kk, ]
      
      if(plot.crit.bic){
        par(mfrow= c(2, 5))
        for(mm in 1:length(crits)){
          plot(list.lambda, beta[, (ncol(x)+mm)], type = "l", lwd=0.7, col='darkblue', xlab = "lambda", log = "x", 
               main = paste0(main, " -- ", crits[mm]));
          points(list.lambda, beta[, (ncol(x)+mm)], type = "p", cex=0.8, pch=16, col='darkblue')
          abline(v=c(cv.fit$lambda.min, cv.fit$lambda.1se), lwd=0.8, col='red')
        }
        
      }
      
      if(method == "cv.lambda.min") {
        optimal = which(fit$lambda==cv.fit$lambda.min)
      }else{
        if(method == "cv.lambda.1se"){
          optimal = which(fit$lambda==cv.fit$lambda.1se)
        }else{
          optimal = which(beta.fitered$HBIC2 == min(beta.fitered$HBIC2))[1]
          #optimal.HBIC3 = which(beta.fitered$HBIC3 == min(beta.fitered$HBIC3))[1]
        }
      }
      
    }else{
      cat("--- overfitting for all lambda ---")
    }
    
    #myCoefs <- coef(fit, s=cv.fit$lambda.min);
    #myCoefs = coef(fit, s=cv.fit$lambda.1se)
    #coefs = as.numeric(myCoefs)[-1]
    #rr.coefs = coefs/max(coefs)
    #coefs[which(rr.coefs<0.05)] = 0
    
    #crit=switch(crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc)
    
    #selected=best.model = which(crit == min(crit))
    
    #ic=c(bic=bic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])
    result = list(coefficients = beta[optimal, c(1:ncol(x))], 
                  optimal = optimal, 
                  coefficients.all = beta)
    
  
  }
  
  return(result)
  
}



###############################
# test the multiple step adaptive elastic-net (msaenet) to select the tuning parameters for both 
# alpha and lambda
# since the package is quite complicated in the design and usage and there is no dramatic improvement (I assume)
# compared the initial paper for adaptive elastic-net 
# At the end this pacakge is abandoned. 
###############################
select.tuning.parameters.using.msaenet = function(x, y)
{
  library("msaenet")
  library(doParallel)
  
  Use.masenet = FALSE
  tune.method = "cv";
  tune.method.nstep = "bic"
  #extract.res.from.msa = TRUE;
  
  alpha.sels = c()
  
  for(n in 1:ncol(y) )
  {
    registerDoParallel(detectCores())
    msa.fit = msaenet(x, y[,n], alphas = seq(0.05, 1.0, 0.05), family='gaussian', 
                      lower.limits = 0, rule = "lambda.min", 
                      tune = tune.method, tune.nsteps = tune.method.nstep,
                      nfolds = 5L, nsteps = 10L,  seed= 1005, parallel = TRUE)
    #plot(msa.fit)
    alpha =  msa.fit$best.alphas[msa.fit$best.step]
    cat("alpha = ", alpha, "\n");
    alpha.sels = c(alpha.sels, alpha)
    myCoefs = as.numeric(msa.fit$beta)
    res[,n] = as.numeric(myCoefs)
  }
  
  #if(extract.res.from.msa){
  
  library("pheatmap")
  library("RColorBrewer")
  cols = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100)
  pheatmap(log2(res + 2^-10), cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=FALSE, main = paste0("gene-specific alphas "),
           color = cols)
}


########################################################
########################################################
# Section: OLS for coarse neuron groups
########################################################
########################################################
run.glmnet.for.coarse.groups = function(x, y, method = "lasso", pval.cutoff=0.1)
{
  res = matrix(0, nrow = ncol(x), ncol = ncol(y)) 
  colnames(res) = colnames(y)
  rownames(res) = colnames(x)
  
  if(method == "ols"){
    cat("fitting the data with ols----")
    
    for(n in 1:ncol(y))
    {
      # n = 59;pval.cutoff=0.5
      y.sel = y[, n]
      x.y.sel = data.frame(x, y.sel)
      model =  eval(parse(text = paste0("y.sel ~ ", paste0(colnames(x), collapse = " + "))))
      
      ff = summary(lm(model, data = x.y.sel))
      coefs = coefficients(ff)[-1, ] # ignore the intercept
      kk = which(coefs[,4]<pval.cutoff)
      res[kk, n] = coefs[kk,1]
    }
  }
  
  if(method == "lasso"){
    res = run.glmnet.select.tuning.parameters(x, y, alphas = c(1), method = "cv.lambda.1se", lambda = 10^seq(-3, 3, length.out = 1000),
                                              plot.cluster.col = TRUE);
    res = res[[1]];
  }
  
  return(res)
  
}

########################################################
########################################################
# Section: test the gdcnet for elastic-net
########################################################
########################################################
########################################################
########################################################
# Section: gcdnet with global lambda2 parameter
########################################################
########################################################
TEST.gcdnet.global.lambda2 = FALSE
if(TEST.gcdnet.global.lambda2){
  require(gcdnet)
  library("pheatmap")
  library("RColorBrewer")
  
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






