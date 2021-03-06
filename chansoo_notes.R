rm(list=ls())
load('data/hw4.rdata')

# Test 
confounders = c('bw','b.marr')
df = hw4[hw4$bw<3000, c('ppvtr.36', 'treat', confounders)]

#########################
ps.m1 <- glm(treat ~ ., data=df[,c('treat',confounders)], family=binomial(link='logit'))
df$psc <- ps.m1$fitted.values
matches <- matching(z=df$treat, score=df$psc, replace=TRUE)
wts <- matches$cnts

#########################

examine_balance = function(df, confounders, wts){
  df2 <- df[, confounders]
  binary = apply(df2,2,function(x) { all(x %in% 0:1) }) # Binary Variable Indicator
  trt = df$treat == 1    # Treatment Indicator
  ctr = df$treat == 0    # Control Indicator
  n_trt = sum(trt)       # Treatment Sample Size
  n_ctrl = sum(ctr)      # Control Sample Size
  wts.trt = rep(1,n_trt) # Set treatment weights
  wts.ctr = wts          # Set control weights
  
  # Means
  trt.means = apply(df2[trt,],2,mean)
  trt.means_w = apply(df2[trt,],2,weighted.mean,wts.trt)
  ctr.means = apply(df2[ctr,],2,mean)
  ctr.means_w = apply(df2[ctr,],2,weighted.mean,wts.ctr)
  
  # Variances
  trt.var = apply(df2[trt,],2,var)
  trt.var_w = apply(df2[trt,],2, function(x) sum(wts.trt * (x - weighted.mean(x, wts.trt))^2) / (sum(wts.trt) - 1))
  ctr.var = apply(df2[ctr,],2,var)
  ctr.var_w = apply(df2[ctr,],2, function(x) sum(wts.ctr * (x - weighted.mean(x, wts.ctr))^2) / (sum(wts.ctr) - 1))
  
  # Standardized Mean Differences
  mean.diff = (trt.means-ctr.means) / sqrt(trt.var)
  mean.diff.bin = (trt.means-ctr.means)
  diff = mean.diff*(1-binary) + mean.diff.bin*binary
  
  mean.diff_w = (trt.means_w-ctr.means_w) / sqrt(trt.var_w)
  mean.diff.bin_w = (trt.means_w-ctr.means_w)
  diff.m = mean.diff_w*(1-binary) + mean.diff.bin_w*binary
  
  # Variance Ratios
  ratio = sqrt(ctr.var)/sqrt(trt.var)
  ratio.m = sqrt(ctr.var_w)/sqrt(trt.var_w)
  
  # Print
  result = cbind(trt.means,ctr.means,trt.means,ctr.means_w,diff,diff.m,ratio,ratio.m)
  colnames(result) = c('mn1','mn0','mn1.m','mn0.m','diff','diff.m','ratio','ratio.m')
  return(round(result,3))
}

examine_balance(df,confounders,wts)
