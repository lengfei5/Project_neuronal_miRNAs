##########################################################################
##########################################################################
# Project:
# Script purpose: to select the optimal model using BIC, AIC, AICc and HQC for elastic-net
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 19 11:09:08 2018
##########################################################################
##########################################################################
feature.selection.and.refitting.ols = function(fit, x, y, mse.overfitting.threshold = 10^-5)
{
  # y = y[,n]
  n.data = length(y)
  #n=length(y)
  list.lambda = fit$lambda;
  crits = c("mse", "BIC", "AIC", "AICc", "HQC","EBIC", "HBIC",  "HBIC2", "HBIC3")
  beta = matrix(0, nrow = length(list.lambda), ncol = (ncol(x)+length(crits)))
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
      beta[n, feature.sels] = as.numeric(coefficients(ff))[-1]
    }
    
    nvar = fit$df[n] + 1
    p_full = ncol(x)
    beta[n, (ncol(x)+1)] = mse
    beta[n, (ncol(x)+2)] = n.data*log(mse) + nvar*log(n.data) # bic
    beta[n, (ncol(x)+3)] = n.data*log(mse) + 2*nvar # aic
    beta[n, (ncol(x)+4)] = n.data*log(mse) + 2*nvar + (2*nvar*(nvar+1))/(n.data-nvar-1) # aicc
    beta[n, (ncol(x)+5)] = n.data*log(mse) + 2*nvar*log(log(n.data))
    gamma = 0.5; beta[n, (ncol(x)+6)] = n.data*log(mse) + nvar*(log(n.data) + 2*gamma*log(p_full)) # ebic
    gamma = 0.5; beta[n, (ncol(x)+7)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic
    gamma = 1; beta[n, (ncol(x)+8)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic2
    gamma = 2; beta[n, (ncol(x)+9)] = n.data*log(mse) + 2*gamma*nvar*log(p_full)# hbic3
  }
  
  beta = data.frame(beta)
  #plot(fit, label = TRUE)
  #plot(fit, xvar = "lambda", label = TRUE); abline(v=log(cv.fit$lambda.min))
  kk = which(beta$mse>mse.overfitting.threshold)
  if(length(kk)>0){
    beta.fitered = beta[kk, ]
    
    par(mfrow= c(2,2))
    plot(beta.fitered$BIC)
    #plot(beta.fitered$EBIC)
    plot(beta.fitered$HBIC)
    plot(beta.fitered$HBIC2)
    plot(beta.fitered$HBIC3)
    
    optimal = which(beta.fitered$HBIC2 == min(beta.fitered$HBIC2))[1]
    #optimal.HBIC3 = which(beta.fitered$HBIC3 == min(beta.fitered$HBIC3))[1]
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
                 coefficients.all = beta )
  
  #result=list(coefficients=coef[,selected],ic=ic,lambda = lambda[selected], nvar=nvar[selected],
  #            glmnet=model,residuals=residuals[,selected],fitted.values=yhat[,selected],ic.range=crit, call = match.call())
  
  #class(result)="ic.glmnet"
  return(result)  
  
}
