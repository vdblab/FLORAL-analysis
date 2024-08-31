### This code is meant to be initiated by a shell script with specified simulation parameters
### n: sample size
### p: number of features
### sparsity: proportion of zero features (between 0 and 1)
### rho: correlation between features (between 0 and 1)
### ns: effect size

options(echo=TRUE) # if you want to see commands in output file
args = commandArgs(trailingOnly=T)
n = as.numeric(args[1])
p = as.numeric(args[2])
sparsity = as.numeric(args[3])
rho = as.numeric(args[4])
ns = as.numeric(args[5])

### As an example, the following set up can be used:
# n = 100
# p = 200
# sparsity = 0.8
# rho = 0
# ns = 0.5

library(FLORAL)
library(glmnet)
library(metagenomeSeq)
library(corncob)
library(phyloseq)
library(ANCOMBC)
library(LDM)
library(lefser)
library(readr)
library(survival)
library(compositions)
library(ALDEx2)
library(zeroSum)
source("https://raw.githubusercontent.com/FrederickHuangLin/ANCOM/8abbfb67ed1cb9ae366ab6349e0ad8f017ca6512/scripts/ancom_v2.1.R")

Nsim=100

TPR <- matrix(NA,nrow=Nsim,ncol=21)
FPR <- matrix(NA,nrow=Nsim,ncol=21)

for (i in 1:Nsim){
  
  set.seed(i)
  
  print(i)
  
  # This example uses survival outcome as an example. For other types of outcomes, change "model" argument and "family" argument in simulation and model fitting functions, correspondingly.
  
  simdat <- simu(n=n,p=p,model="cox",pct.sparsity = sparsity,rho=rho,weaksize = 0.5*ns,strongsize = ns)
  countmat <- exp(simdat$x)-1
  relabmat <- apply(countmat,2,function(x) x/rowSums(countmat))
  clrmat <- matrix(clr(exp(simdat$x)),nrow=n,ncol=p)
  yd <- as.numeric(simdat$d==1)
  
  ##### standard un-constrained lasso with log scaled covariates
  
  glmnetfit <- cv.glmnet(x=simdat$x,y=Surv(simdat$t,simdat$d),family="cox",type.measure = "deviance",nfolds=10,keep = TRUE)
  
  glmnetidx.min <- which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[1]]!=0)
  glmnetidx.1se <- which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[2]]!=0)
  TPR[i,1] <- mean(simdat$idx %in% glmnetidx.min)
  FPR[i,1] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.min)
  TPR[i,2] <- mean(simdat$idx %in% glmnetidx.1se)
  FPR[i,2] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.1se)
  
  ##### standard un-constrained lasso with relative abundance covariates
  
  glmnetfit2 <- cv.glmnet(x=relabmat,Surv(simdat$t,simdat$d),family="cox",type.measure = "deviance",foldid=glmnetfit$foldid)
  
  glmnetidx.min <- which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[1]]!=0)
  glmnetidx.1se <- which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[2]]!=0)
  TPR[i,3] <- mean(simdat$idx %in% glmnetidx.min)
  FPR[i,3] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.min)
  TPR[i,4] <- mean(simdat$idx %in% glmnetidx.1se)
  FPR[i,4] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.1se)
  
  # standard un-constrained lasso with clr transformed covariates
  
  glmnetfit3 <- cv.glmnet(x=clrmat,y=Surv(simdat$t,simdat$d),family="cox",type.measure = "deviance",foldid=glmnetfit$foldid)
  
  glmnetidx.min <- which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[1]]!=0)
  glmnetidx.1se <- which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[2]]!=0)
  TPR[i,5] <- mean(simdat$idx %in% glmnetidx.min)
  FPR[i,5] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.min)
  TPR[i,6] <- mean(simdat$idx %in% glmnetidx.1se)
  FPR[i,6] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.1se)
  
  ##### FLORAL
  
  LRRfit <- try(FLORAL(x=simdat$xcount,y=Surv(simdat$t,simdat$d),family="cox",ncv=10,progress=FALSE,plot=FALSE,foldid=glmnetfit$foldid),silent = TRUE)
  if (!("try-error" %in% class(LRRfit))){
    
    if (!is.null(LRRfit$selected.feature$min.2stage)){
      LRRidx.min.step2 <- parse_number(unique(as.vector(LRRfit$selected.feature$min.2stage)))
      TPR[i,7] <- mean(simdat$idx %in% LRRidx.min.step2)
      FPR[i,7] <- mean(setdiff(1:p,simdat$idx) %in% LRRidx.min.step2)
    }
    
    if (!is.null(LRRfit$selected.feature$`1se.2stage`)){
      LRRidx.1se.step2 <- parse_number(unique(as.vector(LRRfit$selected.feature$`1se.2stage`)))
      TPR[i,8] <- mean(simdat$idx %in% LRRidx.1se.step2)
      FPR[i,8] <- mean(setdiff(1:p,simdat$idx) %in% LRRidx.1se.step2)
    }
  }
  
  ##### zeroSum
  
  zerosumfit <- zeroSum(x=simdat$x,y=Surv(simdat$t,simdat$d),family="cox",nfolds=10,foldid=glmnetfit$foldid)
  
  glmnetidx.min <- which(as.matrix(zerosumfit$coef[[zerosumfit$lambdaMinIndex]]) != 0)-1
  glmnetidx.1se <- which(as.matrix(zerosumfit$coef[[zerosumfit$lambda1SEIndex]]) != 0)-1
  TPR[i,9] <- mean(simdat$idx %in% glmnetidx.min)
  FPR[i,9] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.min)
  TPR[i,10] <- mean(simdat$idx %in% glmnetidx.1se)
  FPR[i,10] <- mean(setdiff(1:p,simdat$idx) %in% glmnetidx.1se)
  
  ##### metagenomeSeq
  
  test_obj<- newMRexperiment(
    counts = t(exp(simdat$x)-1), 
    phenoData = AnnotatedDataFrame(data.frame(outcome=yd))
  )
  pp <- cumNormStat(test_obj, pFlag = T) #Cumulative sum scaling percentile selection
  test_obj_norm <- cumNorm(test_obj, p=pp) #Cumulative sum scaling normalization
  mod_class <- model.matrix(
    ~outcome, 
    data=pData(test_obj_norm)
  )
  regres_class <-  try(fitZig(test_obj_norm, mod_class),silent=TRUE)
  if ("try-error" %in% class(regres_class)){
    TPR[i,11] <- NA
    FPR[i,11] <- NA
  }else{
    res_table_class <- MRfulltable(regres_class, coef="outcome", number = p)
    metagenomeSeq_hits <- rownames(res_table_class)[res_table_class$adjPvalues < 0.05]
    TPR[i,11] <- mean(simdat$idx %in% metagenomeSeq_hits)
    FPR[i,11] <- mean(setdiff(1:p,simdat$idx) %in% metagenomeSeq_hits)
  }
  
  ##### corncob
  
  OTU = otu_table(exp(simdat$x)-1, taxa_are_rows = FALSE)
  META = sample_data(data.frame(outcome=as.factor(yd)))
  physeq = phyloseq(OTU, META)
  da_analysis <- try(differentialTest(
    formula =  ~outcome,
    formula_null = ~1,
    phi.formula = ~1,
    phi.formula_null = ~1,
    test = "Wald", boot = FALSE,
    data = physeq,
    fdr_cutoff = 0.05, verbose=TRUE,
    control = list(maxit = 1000, reltol = 1e-1)
  ),silent = TRUE)
  if (!("try-error" %in% class(da_analysis))){
    TPR[i,12] <- mean(simdat$idx %in% da_analysis$significant_taxa)
    FPR[i,12] <- mean(setdiff(1:p,simdat$idx) %in% da_analysis$significant_taxa)
  }
  
  ##### ANCOM-BC
  
  OTU = otu_table(exp(simdat$x)-1, taxa_are_rows = FALSE)
  META = sample_data(data.frame(outcome=as.factor(yd)))
  physeq = phyloseq(OTU, META)
  res = try(ancombc(
    phyloseq = physeq, 
    formula = "outcome",
    p_adj_method = "BH", 
    lib_cut = 0, 
    group = "outcome", 
    struc_zero = TRUE, 
    neg_lb = TRUE, 
    tol = 1e-5, 
    max_iter = 100, 
    conserve = TRUE, 
    alpha = 0.05,
    global = TRUE),silent=TRUE)
  
  if (!("try-error" %in% class(res))){
    TPR[i,13] <- mean(simdat$idx %in% which(res$res$diff_abn == "TRUE"))
    FPR[i,13] <- mean(setdiff(1:p,simdat$idx) %in% which(res$res$diff_abn == "TRUE"))
  }
  
  ##### LDM
  
  cox.fit = coxph(Surv(simdat$t, simdat$d==1) ~ 1)
  resid.Martingale = residuals(cox.fit, type='martingale')
  clin <- data.frame(outcome=resid.Martingale)
  otudat <- exp(simdat$x)-1
  res.ldm <- try(ldm(otudat ~ resid.Martingale, data=clin, seed=41222),silent=TRUE)
  
  if ("try-error" %in% class(res.ldm)){
    TPR[i,14] <- NA
    FPR[i,14] <- NA
  }else{
    TPR[i,14] <- mean(simdat$idx %in% as.numeric(res.ldm$detected.otu.omni[[1]]))
    FPR[i,14] <- mean(setdiff(1:p,simdat$idx) %in% as.numeric(res.ldm$detected.otu.omni[[1]]))
  }
  
  ##### LEFSE
  
  se0 <- SummarizedExperiment(assays=SimpleList(counts=t(exp(simdat$x)-1)),
                              colData=data.frame(outcome=as.factor(yd)))
  res_group <- try(lefser(se0, groupCol = "outcome",
                          kruskal.threshold = 0.05,
                          wilcox.threshold = 0.05),silent = TRUE)
  
  if (!("try-error" %in% class(res_group))){
    TPR[i,15] <- mean(simdat$idx %in% parse_number(res_group$Names))
    FPR[i,15] <- mean(setdiff(1:p,simdat$idx) %in% parse_number(res_group$Names))
  }
  
  ##### ANCOM-II
  
  meta_data <- data.frame(outcome=yd,id=1:n)
  rownames(simdat$x) <- 1:n
  prepro = feature_table_pre_process(
    feature_table = t(exp(simdat$x)-1), meta_data=meta_data, sample_var="id",
    out_cut= 0.05, zero_cut= 0.90, lib_cut = 0, neg_lb=TRUE)
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  rownames(prepro$feature_table) <- rownames(prepro$feature_table)
  
  res = try(ANCOM(
    feature_table = prepro$feature_table,
    meta_data = prepro$meta_data,
    struc_zero=prepro$structure_zeros ,# Structural zero info,
    main_var = "outcome", 
    p_adj_method="BH", 
    alpha = 0.05,
    adj_formula=NULL,
    rand_formula=NULL, 
    control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")),silent=TRUE)
  
  if (!("try-error" %in% class(res))){
    TPR[i,16] <- mean(simdat$idx %in% as.numeric(res$out$taxa_id[res$out$detected_0.9]))
    FPR[i,16] <- mean(setdiff(1:p,simdat$idx) %in% as.numeric(res$out$taxa_id[res$out$detected_0.9]))
  }
  
  ##### Linear regression / Wilcoxon test with taxa as outcome
  
  ######### Relative abundance
  
  lmmat <- matrix(NA,nrow=p,ncol=2)
  wimat <- matrix(NA,nrow=p,ncol=2)
  
  for (j in 1:p){
    fit <- summary(lm(log(relabmat[,j]+1e-8)~yd))
    lmmat[j,1] <- fit$coefficients[2,4]
    wil <- wilcox.test(relabmat[yd==1,j],relabmat[yd==0,j])
    wimat[j,1] <- wil$p.value
  }
  
  lmmat[,2] <- p.adjust(lmmat[,1])
  TPR[i,17] <- mean(simdat$idx %in% which(lmmat[,2]<0.05))
  FPR[i,17] <- mean(setdiff(1:p,simdat$idx) %in% which(lmmat[,2]<0.05))
  
  wimat[,2] <- p.adjust(wimat[,1])
  TPR[i,18] <- mean(simdat$idx %in% which(wimat[,2]<0.05))
  FPR[i,18] <- mean(setdiff(1:p,simdat$idx) %in% which(wimat[,2]<0.05))
  
  ######## Centered log-ratio
  
  lmmat <- matrix(NA,nrow=p,ncol=2)
  wimat <- matrix(NA,nrow=p,ncol=2)
  
  for (j in 1:p){
    fit <- summary(lm(clrmat[,j]~yd))
    lmmat[j,1] <- fit$coefficients[2,4]
    wil <- wilcox.test(clrmat[yd==1,j],clrmat[yd==0,j])
    wimat[j,1] <- wil$p.value
  }
  lmmat[,2] <- p.adjust(lmmat[,1])
  TPR[i,19] <- mean(simdat$idx %in% which(lmmat[,2]<0.05))
  FPR[i,19] <- mean(setdiff(1:p,simdat$idx) %in% which(lmmat[,2]<0.05))
  
  wimat[,2] <- p.adjust(wimat[,1])
  TPR[i,20] <- mean(simdat$idx %in% which(wimat[,2]<0.05))
  FPR[i,20] <- mean(setdiff(1:p,simdat$idx) %in% which(wimat[,2]<0.05))
  
  ##### ALDEx2
  
  conds <- yd
  mat <- t(exp(simdat$x)-1)
  mode(mat) <- "integer"
  res = try(aldex(mat, 
                  conds, 
                  mc.samples=1000, 
                  test="t", 
                  effect=TRUE,
                  include.sample.summary=FALSE, 
                  denom="all", 
                  verbose=FALSE,
                  paired.test=FALSE)
            ,silent=TRUE)
  if (!("try-error" %in% class(res))){
    TPR[i,21] <- mean(simdat$idx %in% which(res$we.eBH < 0.05))
    FPR[i,21] <- mean(setdiff(1:p,simdat$idx) %in% which(res$we.eBH < 0.05))
  }
  
  gc()
  
}

F1score <- 2*TPR*10/(2*TPR*10 + FPR*(p-10) + 10-TPR*10)
FP <- FPR*(p-10)
FN <- 10-TPR*10

colnames(F1score) <- colnames(FP) <- colnames(FN) <-  
  c("glmnet(log), lambda(min)",
    "glmnet(log), lambda(1se)",
    "glmnet(relative abundance), lambda(min)",
    "glmnet(relative abundance), lambda(1se)",
    "glmnet(clr), lambda(min)",
    "glmnet(clr), lambda(1se)",
    "FLORAL, lambda(min)",
    "FLORAL, lambda(1se)",
    "zeroSum, lambda(min)",
    "zeroSum, lambda(1se)",
    "metagenomeSeq",
    "corncob",
    "ANCOM-BC",
    "LDM",
    "LEfSe",
    "ANCOM-II",
    "LM(relative abundance)",
    "Wilcoxon(relative abundance)",
    "LM(clr)",
    "Wilcoxon(clr)",
    "ALDEx2")

