feature.selection.and.refitting.ols = function(fit, x, y)
{
  # y = y[,n]
  nb.data = length(y)
  #n=length(y)
  list.lambda = fit$lambda;
  bic=rep(NA, length(list.lambda));
  aic=bic;
  aicc = bic; 
  hqc=bic;
  for(n in 1:length(list.lambda))
  {
    # n = 20;
    #df = fit$df[n]
    cfs = as.numeric(coef(fit, s=list.lambda[n]))
    cfs = cfs[-1] # move the intercept
    feature.sels = which(cfs>0)
    names(feature.sels) = colnames(x)[feature.sels]
    if(length(feature.sels)==0) {
      mse = sum(y-mean(y)^2);
    }else{
      x.sel = x[, feature.sels]
      ff = lm(y ~ x.sel)
      mse = sum(ff$residuals^2)
    }
    
    nvar = fit$df[n] + 1
    bic = nb.data*log(mse)+nvar*log(nb.data)
    aic = nb.data*log(mse)+2*nvar
    aicc = aic+(2*nvar*(nvar+1))/(nb.data-nvar-1)
    hqc = nb.data*log(mse)+2*nvar*log(log(nb.data))
    
  }
  
  #plot(fit, label = TRUE)
  #plot(fit, xvar = "lambda", label = TRUE); abline(v=log(cv.fit$lambda.min))
  
  #myCoefs <- coef(fit, s=cv.fit$lambda.min);
  #myCoefs = coef(fit, s=cv.fit$lambda.1se)
  #coefs = as.numeric(myCoefs)[-1]
  #rr.coefs = coefs/max(coefs)
  #coefs[which(rr.coefs<0.05)] = 0
  
  #model = glmnet(x = x, y = y, ...)
  #coef = coef(model)
  #lambda = model$lambda
  #df = model$df
  
  #yhat=cbind(1,x)%*%coef
  #residuals = (y- yhat)
  #mse = colMeans(residuals^2)
  #sse = colSums(residuals^2)
  
 
  
  sst = (n-1)*var(y)
  r2 = 1 - (sse/sst)
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))
  
  crit=switch(crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc)
  
  selected=best.model = which(crit == min(crit))
  
  ic=c(bic=bic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])
  
  result=list(coefficients=coef[,selected],ic=ic,lambda = lambda[selected], nvar=nvar[selected],
              glmnet=model,residuals=residuals[,selected],fitted.values=yhat[,selected],ic.range=crit, call = match.call())
  
  class(result)="ic.glmnet"
  return(result)  
  
}
