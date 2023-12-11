### Please download the dataset folder "Hackathon" from https://figshare.com/articles/dataset/16S_rRNA_Microbiome_Datasets/14531724

### Set the path to the folder Hackathon
datapath="Hackathon"

### Before running this code, please create a folder named "results_seed1" under the folder "Hackathon" to save results.

library(FLORAL)
library(tidyverse)
library(glmnet)
library(pROC)
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

studies = list.dirs(file.path(datapath, "Studies"), recursive = FALSE)

set.seed(1) # Set seed of shuffling y

# The following loop will run for all datasets included. 
# Please specify the index of the data set to run specific ones using code inside the loop.

for (i in 1:length(studies)){ 
  
  study = studies[i]
  print(study)
  
  ### Load metadata table
  metadata <- read.csv(sep="\t",list.files(study,pattern = "_meta", full.names = TRUE), check.names = FALSE, skip=1, header=FALSE) %>% dplyr::rename(SampleID=1,y=2) %>% mutate(SampleID = as.character(SampleID))
  
  
  ### Load genus table
  data_path <- list.files(study,pattern = "_genus_table.tsv", full.names = TRUE)
  # some data files have a extra header so we read the header in separate from the rest of the data
  these_lines <- readLines(data_path)
  sample_list <- strsplit(rev(these_lines[grepl("^#", these_lines)])[1], "\t")[[1]]
  sample_list <- sample_list[sample_list != "#OTU ID"]
  genusdata <- read.csv(sep="\t", data_path, comment.char = "#",header=FALSE) %>% column_to_rownames("V1") %>%  t() 
  
  rownames(genusdata) <- sample_list
  both <- inner_join(metadata, genusdata %>% data.frame() %>% rownames_to_column("SampleID")) 
  
  ### Derive the binary outcome as y, and feature count matrix as x
  
  y = as.numeric(as.factor(both$y)) - 1
  y = sample(y,size=length(y)) ### Shuffle the label y randomly, for results without shuffling, comment this line
  x = both %>% dplyr::select(-y, -SampleID) %>% as.matrix()
  
  Selected <- data.frame()
  
  ### glmnet log(x)
  
  tstart <- proc.time()
  glmnetfit <- cv.glmnet(x=log1p(x),y=y,family="binomial",type.measure = "mse",nfolds=5,keep=TRUE)
  tstop <- proc.time()
  
  if (length(which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[1]]!=0)) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=colnames(x)[which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[1]]!=0)],
                                 method="glmnet, log, min",
                                 coef=glmnetfit$glmnet.fit$beta[,glmnetfit$index[1]][which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[1]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, log, min",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  if (length(which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[2]]!=0))){
    Selected <- rbind(Selected, 
                      data.frame(taxon=colnames(x)[which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[2]]!=0)],
                                 method="glmnet, log, 1se",
                                 coef=glmnetfit$glmnet.fit$beta[,glmnetfit$index[2]][which(glmnetfit$glmnet.fit$beta[,glmnetfit$index[2]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, log, 1se",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### standard un-constrained lasso with clr transformed covariates
  
  tstart <- proc.time()
  glmnetfit2 <- cv.glmnet(x=clrmat,y=y,family="binomial",type.measure = "mse",nfolds=5,foldid = glmnetfit$foldid)
  tstop <- proc.time()

  if (length(which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[1]]!=0)) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon= colnames(x)[which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[1]]!=0)],
                                 method="glmnet, clr, min",
                                 coef=glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[1]][which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[1]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, clr, min",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  if (length(which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[2]]!=0)) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon= colnames(x)[which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[2]]!=0)],
                                 method="glmnet, clr, 1se",
                                 coef=glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[2]][which(glmnetfit2$glmnet.fit$beta[,glmnetfit2$index[2]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, clr, 1se",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### standard un-constrained lasso with relative abundance covariates
  
  tstart <- proc.time()
  glmnetfit3 <- cv.glmnet(x=relabmat,y=y,family="binomial",type.measure = "mse",nfolds=5,foldid = glmnetfit$foldid)
  tstop <- proc.time()
  
  if (length(which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[1]]!=0)) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon= colnames(x)[which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[1]]!=0)],
                                 method="glmnet, rel, min",
                                 coef=glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[1]][which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[1]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, rel, min",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  if (length(which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[2]]!=0)) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon= colnames(x)[which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[2]]!=0)],
                                 method="glmnet, rel, 1se",
                                 coef=glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[2]][which(glmnetfit3$glmnet.fit$beta[,glmnetfit3$index[2]]!=0)],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="glmnet, rel, 1se",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### FLORAL
  
  tstart <- proc.time()
  FLORALfit <- try(FLORAL(x,y,family="binomial",a=1,ncv=5,step2=TRUE,plot=FALSE,progress = TRUE,foldid=glmnetfit$foldid))
  tstop <- proc.time()
  
  if (class(FLORALfit) != "try-error"){
    
    if (length(FLORALfit$selected.feature$min.2stage) > 0){
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=FLORALfit$selected.feature$min.2stage,
                                   method="FLORAL, min 2 stage",
                                   coef=FLORALfit$best.beta$min[FLORALfit$selected.feature$min.2stage],
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }else{
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=NA,
                                   method="FLORAL, min 2 stage",
                                   coef=NA,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }
    
    if (length(FLORALfit$selected.feature$`1se.2stage`) > 0){
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=FLORALfit$selected.feature$`1se.2stage`,
                                   method="FLORAL, 1se 2 stage",
                                   coef=FLORALfit$best.beta$`1se`[FLORALfit$selected.feature$`1se.2stage`],
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
    }else{
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=NA,
                                   method="FLORAL, 1se 2 stage",
                                   coef=NA,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }
    
  }
  
  ### zeroSum
  
  tstart <- proc.time()
  zerosumfit <- zeroSum(x=log1p(x),y=y,family="binomial",nfolds=5,foldid = glmnetfit$foldid)
  tstop <- proc.time()

  if (length(which(as.matrix(zerosumfit$coef[[zerosumfit$lambdaMinIndex]]) != 0)) > 1){
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=colnames(x)[which(as.matrix(zerosumfit$coef[[zerosumfit$lambdaMinIndex]]) != 0)-1],
                                 method="zeroSum, min",
                                 coef=as.matrix(zerosumfit$coef)[[zerosumfit$lambdaMinIndex]][which(as.matrix(zerosumfit$coef[[zerosumfit$lambdaMinIndex]]) != 0)][-1],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="zeroSum, min",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  if (length(which(as.matrix(zerosumfit$coef[[zerosumfit$lambda1SEIndex]]) != 0)) > 1){
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=colnames(x)[which(as.matrix(zerosumfit$coef[[zerosumfit$lambda1SEIndex]]) != 0)-1],
                                 method="zeroSum, 1se",
                                 coef=as.matrix(zerosumfit$coef)[[zerosumfit$lambda1SEIndex]][which(as.matrix(zerosumfit$coef[[zerosumfit$lambda1SEIndex]]) != 0)][-1],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="zeroSum, 1se",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### LDM
  
  clin <- data.frame(outcome=y)
  otudat <- x
  tstart <- proc.time()
  res.ldm <- try(ldm(otudat ~ outcome, 
                     data=clin, 
                     fdr.nominal = 0.05,
                     seed=41222),silent=TRUE)
  tstop <- proc.time()

  if (class(res.ldm) != "try-error"){
    
    if (length(res.ldm$detected.otu.omni[[1]]) > 0){
      Selected <- rbind(Selected, 
                        data.frame(taxon=res.ldm$detected.otu.omni[[1]],
                                   method="LDM",
                                   coef=as.vector(t(res.ldm$beta[,res.ldm$detected.otu.omni[[1]]])),
                                   pval=res.ldm$q.otu.omni[,res.ldm$detected.otu.omni[[1]]],
                                   time=tstop[3]-tstart[3])
      )
    }else{
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=NA,
                                   method="LDM",
                                   coef=NA,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }
  }
  
  relabmat <- apply(x,2,function(y) y/rowSums(x))
  clrmat <- matrix(clr(x),nrow=nrow(x),ncol=ncol(x))
  
  
  ### LEFSE

  se0 <- SummarizedExperiment(assays=SimpleList(counts=t(x)),
                              colData=data.frame(outcome=as.factor(y)))
  
  tstart <- proc.time()
  res_group <- try(lefser(se0, groupCol = "outcome",
                          kruskal.threshold = 0.05,
                          wilcox.threshold = 0.05),silent = TRUE)
  tstop <- proc.time()

  if (class(res_group) != "try-error"){
    if (length(res_group$Names) > 0){
      Selected <- rbind(Selected, 
                        data.frame(taxon=res_group$Names,
                                   method="LEfSe",
                                   coef=res_group$scores,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
    }
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="LEfSe",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### ANCOM-II
  
  meta_data <- data.frame(outcome=y,id=1:length(y))
  rownames(x) <- 1:length(y)
  
  tstart <- proc.time()
  prepro = feature_table_pre_process(
    feature_table = t(x), meta_data=meta_data, sample_var="id",
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
  tstop <- proc.time()
  
  if (length(res$out$taxa_id[res$out$detected_0.9]) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=res$out$taxa_id[res$out$detected_0.9],
                                 method="ANCOM",
                                 coef=res$fig$data[res$out$taxa_id[res$out$detected_0.9],2],
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="ANCOM",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### Linear regression / Wilcoxon test with taxa as outcome (relative abundance)
  
  lmmat <- matrix(NA,nrow=ncol(x),ncol=3)
  wimat <- matrix(NA,nrow=ncol(x),ncol=2)
  rownames(lmmat) <- rownames(wimat) <- colnames(x)
  
  tstart <- proc.time()
  for (j in 1:ncol(x)){
    fit <- summary(lm(log(relabmat[,j]+1e-8)~y))
    lmmat[j,1] <- fit$coefficients[2,4]
    lmmat[j,3] <- coef(fit)[2,1]
  }
  lmmat[,2] <- p.adjust(lmmat[,1])
  tstop <- proc.time()
  
  if (length(rownames(lmmat)[lmmat[,2] < 0.05]) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=rownames(lmmat)[lmmat[,2] < 0.05],
                                 method="lm",
                                 coef=lmmat[lmmat[,2] < 0.05,3],
                                 pval=lmmat[lmmat[,2] < 0.05,2],
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="lm",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  tstart <- proc.time()
  for (j in 1:ncol(x)){
    wil <- wilcox.test(relabmat[y==1,j],relabmat[y==0,j])
    wimat[j,1] <- wil$p.value
  }
  wimat[,2] <- p.adjust(wimat[,1])
  tstop <- proc.time()
  
  if (length(rownames(wimat)[wimat[,2] < 0.05]) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=rownames(wimat)[wimat[,2] < 0.05],
                                 method="Wilcoxon",
                                 coef=NA,
                                 pval=wimat[wimat[,2] < 0.05,2],
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="Wilcoxon",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### ALDEx
  
  if (FALSE){
    
    conds <- y
    mat <- t(x)
    mode(mat) <- "integer"
    
    tstart <- proc.time()
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
    tstop <- proc.time()

    if (length(which(res$we.eBH < 0.05)) > 0){
      Selected <- rbind(Selected, 
                        data.frame(taxon=colnames(x)[which(res$we.eBH < 0.05)],
                                   method="ALDEx",
                                   coef=res$effect[which(res$we.eBH < 0.05)],
                                   pval=res$we.eBH[which(res$we.eBH < 0.05)],
                                   time=tstop[3]-tstart[3])
      )
    }else{
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=NA,
                                   method="ALDEx",
                                   coef=NA,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }
    
  }
  
  ### metagenomeSeq
  
  tstart <- proc.time()
  test_obj<- newMRexperiment(
    counts = t(x), 
    phenoData = AnnotatedDataFrame(data.frame(outcome=y))
  )
  pp <- cumNormStat(test_obj, pFlag = T) #Cumulative sum scaling percentile selection
  test_obj_norm <- cumNorm(test_obj, p=pp) #Cumulative sum scaling normalization
  mod_class <- model.matrix(
    ~outcome, 
    data=pData(test_obj_norm)
  )
  regres_class <-  try(fitZig(test_obj_norm, mod_class),silent=TRUE)
  
  if (class(regres_class) == "try-error"){
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="metagenomeSeq",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }else{
    res_table_class <- MRfulltable(regres_class, coef="outcome", number = ncol(x))
    metagenomeSeq_hits <- rownames(res_table_class)[res_table_class$adjPvalues < 0.05]
    tstop <- proc.time()

    if (length(metagenomeSeq_hits) > 0){
      Selected <- rbind(Selected, 
                        data.frame(taxon=metagenomeSeq_hits,
                                   method="metagenomeSeq",
                                   coef=log(res_table_class$oddsRatio[res_table_class$adjPvalues < 0.05]+0.001),
                                   pval=res_table_class$adjPvalues[res_table_class$adjPvalues < 0.05],
                                   time=tstop[3]-tstart[3])
      )
    }else{
      
      Selected <- rbind(Selected, 
                        data.frame(taxon=NA,
                                   method="metagenomeSeq",
                                   coef=NA,
                                   pval=NA,
                                   time=tstop[3]-tstart[3])
      )
      
    }
  }
  
  ### corncob
  rownames(x) <- NULL
  OTU = otu_table(x, taxa_are_rows = FALSE)
  META = sample_data(data.frame(outcome=as.factor(y)))
  physeq = phyloseq(OTU, META)

  tstart <- proc.time()
  da_analysis <- try(differentialTest(
    formula =  ~outcome,
    formula_null = ~1,
    phi.formula = ~1,
    phi.formula_null = ~1,
    test = "Wald", boot = FALSE,
    data = physeq,
    fdr_cutoff = 0.05, verbose=FALSE,
    control = list(maxit = 1000, reltol = 1e-1)
  ),silent = TRUE)
  tstop <- proc.time()

  pda <- plot(da_analysis)
  rownames(pda$data) = pda$data$taxa
  
  if (length(da_analysis$significant_taxa) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=da_analysis$significant_taxa,
                                 method="corncob",
                                 coef=pda$data[da_analysis$significant_taxa,1],
                                 pval=da_analysis$p_fdr[da_analysis$significant_taxa],
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="corncob",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  ### ANCOM-BC
  OTU = otu_table(x, taxa_are_rows = FALSE)
  META = sample_data(data.frame(outcome=as.factor(y)))
  physeq = phyloseq(OTU, META)
  
  tstart <- proc.time()
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
  tstop <- proc.time()
  
  if (length(which(res$res$diff_abn$outcome1 == "TRUE")) > 0){
    Selected <- rbind(Selected, 
                      data.frame(taxon=rownames(res$feature_table)[which(res$res$diff_abn$outcome1 == "TRUE")],
                                 method="ANCOMBC",
                                 coef=res$res$beta[rownames(res$feature_table)[which(res$res$diff_abn$outcome1 == "TRUE")],],
                                 pval=res$res$q_val[rownames(res$feature_table)[which(res$res$diff_abn$outcome1 == "TRUE")],],
                                 time=tstop[3]-tstart[3])
    )
  }else{
    
    Selected <- rbind(Selected, 
                      data.frame(taxon=NA,
                                 method="ANCOMBC",
                                 coef=NA,
                                 pval=NA,
                                 time=tstop[3]-tstart[3])
    )
    
  }
  
  Selected$study <- basename(study)
  Selected$n <- length(y)
  Selected$p <- ncol(x)
  
  save(Selected,file = paste0("Hackathon/results_seed1/",basename(study),"seed1schuffle.rdata"))
  
  gc()
}
