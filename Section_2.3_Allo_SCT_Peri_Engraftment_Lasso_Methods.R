options(echo=TRUE)
args = commandArgs(trailingOnly=T)
library(dplyr)
library(tidyverse)
library(lubridate)
library(phyloseq)
library(FLORAL)
library(survival)
library(gtsummary)
library(glmnet)
library(zeroSum)
library(compositions)

### Multiple runs of cross validations are performed by parallel computing in this script
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)

# Load microbiome data collected at peri-engraftment period
load("peri_engraft.Rdata")

# Filter out ASVs not identified from the samples
phy_abs = filter_taxa(phy_peri_engraft_filtered, function(x) mean(x > 0) > 0, TRUE) 

# Aggregate ASV counts at genus level
df_taxa_abs <- left_join(
  tax_table(phy_abs) %>% as.data.frame()  %>% rownames_to_column("key"),
  t(otu_table(phy_abs)) %>% as.data.frame() %>% rownames_to_column("key")
) %>% 
  pivot_longer(cols = -c(key:species)) %>% 
  dplyr::select(-c(key, kingdom, phylum,class,    ordr,  family,   species)) %>% 
  pivot_wider(names_from="genus", values_fn = sum, values_fill = 0) %>% 
  column_to_rownames("name")

# Extract metadata from the phyloseq object
df_meta <- sample_data(phy_abs) %>% data.frame()
df_meta <- df_meta %>% arrange(mrn,day_relative_to_hct) 

# Sort the microbiome count table in the same order as the metadata
df_taxa_abs <- df_taxa_abs[rownames(df_meta),]
sum(rownames(df_taxa_abs)!=rownames(df_meta))

# Collapse covariate levels
df_meta <- df_meta %>% mutate(disease.simple = factor(ifelse(disease %in% c("Leukemia","Leukemia (Second Primary)","Myelodysplastic Syndrome","Myelodysplastic Syndromes (Second Primary)","Myeloproliferative Disorder","Myeloproliferative Disorder (Second Primary)"),"AML/ALL/MDS/MPN","Others"),levels=c("AML/ALL/MDS/MPN","Others")),
                              source.simple = factor(ifelse(source %in% c("Allo BM Unmodified",
                                                                          "Allo PBSC Unmodified Only"),
                                                            "Unmodified",
                                                            ifelse(source %in% c("Allo Double Cord Blood",
                                                                                 "Allo PBSC CD34+/Double Cord Blood",
                                                                                 "Allo PBSC CD34+/SINGLE Cord Blood"),
                                                                   "Cord","TCD")),
                                                     levels=c("TCD","Unmodified","Cord")))

# Break the ties in survival times
df_meta$OStime_landmark <- df_meta$OStime_landmark + 1e-8

# Reformat the taxa table as matrix object
df_taxa_abs <- as.matrix(df_taxa_abs)

# Extract the covariate data frame to be adjusted for
covmat <- df_meta %>% select(disease.simple,
                             source.simple,
                             age,
                             intensity) %>% 
  as.data.frame()

# Reformat the covariate table with dummy variables
covmat <- model.matrix(~disease.simple+source.simple+age+intensity,covmat)[,-1]
colnames(covmat)[6] = "intensityReduced.Intensity"

# Create relative abundance and CLR-transformed taxa matrices
relabmat <- apply(df_taxa_abs,2,function(y) y/rowSums(df_taxa_abs))
clrmat <- matrix(clr(df_taxa_abs),nrow=nrow(df_taxa_abs),ncol=ncol(df_taxa_abs))
colnames(clrmat) <- colnames(df_taxa_abs)

### Run FLORAL with 100 random cross validations, using 10 threads

mcv = 100
ncore= 10
cl <- makeCluster(ncore)
registerDoParallel(cl)
seed <- mcv
registerDoRNG(seed=seed)

FLORAL.res <- foreach(i=1:mcv,.packages = c("FLORAL","survival")) %dopar% {
  
  FLORAL(x=cbind(covmat,df_taxa_abs),
         y=Surv(df_meta$OStime_landmark, df_meta$OSevent),
         family="cox",
         ncov=6,
         a=1,
         progress=FALSE,
         plot=FALSE)
  
}

stopCluster(cl)

### Summarize the taxa selection results

res.FLORAL <- list(min=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min)))/mcv,
                   `1se`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se`)))/mcv,
                   min.2stage=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$min.2stage)))/mcv,
                   `1se.2stage`=table(unlist(lapply(FLORAL.res,function(x) x$selected.feature$`1se.2stage`)))/mcv,
                   min.2stage.ratios=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$min))))/mcv,
                   `1se.2stage.ratios`=table(unlist(lapply(FLORAL.res,function(x) rownames(x$step2.tables$`1se`))))/mcv,
                   mcv=mcv,
                   seed=seed)

res.FLORAL$min.coef = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$min))[names(res.FLORAL$min)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$min), function(x) x != 0))[names(res.FLORAL$min)]
res.FLORAL$`1se.coef` = Reduce(`+`,lapply(FLORAL.res,function(x) x$best.beta$`1se`))[names(res.FLORAL$`1se`)]/Reduce(`+`,lapply(lapply(FLORAL.res, function(x) x$best.beta$`1se`), function(x) x != 0))[names(res.FLORAL$`1se`)]

### Run glmnet (log-transformed counts) with 100 cross validations, using 10 threads
### Note that the fold splits are the same as used in FLORAL for each run

cl <- makeCluster(ncore)
registerDoParallel(cl)

glmnetlog.res <- foreach(i=1:mcv,.packages = c("glmnet","survival")) %dopar% {
  
  cv.glmnet(x=cbind(covmat,log1p(df_taxa_abs)),
            y=Surv(df_meta$OStime_landmark, df_meta$OSevent),
            family="cox",
            penalty.factor=c(rep(0,6),rep(1,105)),
            nfolds=5,
            foldid = FLORAL.res[[i]]$foldid
  )
  
}

stopCluster(cl)

res.glmnetlog <- list(min=table(unlist(lapply(glmnetlog.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[1]]!=0))])))/mcv,
                      `1se`=table(unlist(lapply(glmnetlog.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[2]]!=0))])))/mcv,
                      mcv=mcv,
                      seed=seed)

res.glmnetlog$min.coef = Reduce(`+`,lapply(glmnetlog.res,function(x) x$glmnet.fit$beta[,x$index[1]]))[names(res.glmnetlog$min)]/Reduce(`+`,lapply(lapply(glmnetlog.res, function(x) x$glmnet.fit$beta[,x$index[1]]), function(x) x != 0))[names(res.glmnetlog$min)]
res.glmnetlog$`1se.coef` = Reduce(`+`,lapply(glmnetlog.res,function(x) x$glmnet.fit$beta[,x$index[2]]))[names(res.glmnetlog$`1se`)]/Reduce(`+`,lapply(lapply(glmnetlog.res, function(x) x$glmnet.fit$beta[,x$index[2]]), function(x) x != 0))[names(res.glmnetlog$`1se`)]

### Run glmnet (relative abundance) with 100 cross validations, using 10 threads
### Note that the fold splits are the same as used in FLORAL for each run

cl <- makeCluster(ncore)
registerDoParallel(cl)

glmnetrel.res <- foreach(i=1:mcv,.packages = c("glmnet","survival")) %dopar% {
  
  cv.glmnet(x=cbind(covmat,relabmat),
            y=Surv(df_meta$OStime_landmark, df_meta$OSevent),
            family="cox",
            penalty.factor=c(rep(0,6),rep(1,105)),
            nfolds=5,
            foldid = FLORAL.res[[i]]$foldid
  )
  
}

stopCluster(cl)

res.glmnetrel <- list(min=table(unlist(lapply(glmnetrel.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[1]]!=0))])))/mcv,
                      `1se`=table(unlist(lapply(glmnetrel.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[2]]!=0))])))/mcv,
                      mcv=mcv,
                      seed=seed)

res.glmnetrel$min.coef = Reduce(`+`,lapply(glmnetrel.res,function(x) x$glmnet.fit$beta[,x$index[1]]))[names(res.glmnetrel$min)]/Reduce(`+`,lapply(lapply(glmnetrel.res, function(x) x$glmnet.fit$beta[,x$index[1]]), function(x) x != 0))[names(res.glmnetrel$min)]
res.glmnetrel$`1se.coef` = Reduce(`+`,lapply(glmnetrel.res,function(x) x$glmnet.fit$beta[,x$index[2]]))[names(res.glmnetrel$`1se`)]/Reduce(`+`,lapply(lapply(glmnetrel.res, function(x) x$glmnet.fit$beta[,x$index[2]]), function(x) x != 0))[names(res.glmnetrel$`1se`)]

### Run glmnet (clr-transformed counts) with 100 cross validations, using 10 threads
### Note that the fold splits are the same as used in FLORAL for each run

cl <- makeCluster(ncore)
registerDoParallel(cl)

glmnetclr.res <- foreach(i=1:mcv,.packages = c("glmnet","survival")) %dopar% {
  
  cv.glmnet(x=cbind(covmat,clrmat),
            y=Surv(df_meta$OStime_landmark, df_meta$OSevent),
            family="cox",
            penalty.factor=c(rep(0,6),rep(1,105)),
            nfolds=5,
            foldid = FLORAL.res[[i]]$foldid
  )
  
}

stopCluster(cl)

res.glmnetclr <- list(min=table(unlist(lapply(glmnetclr.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[1]]!=0))])))/mcv,
                      `1se`=table(unlist(lapply(glmnetclr.res,function(x) colnames(cbind(covmat,log1p(df_taxa_abs)))[which(as.matrix(x$glmnet.fit$beta[,x$index[2]]!=0))])))/mcv,
                      mcv=mcv,
                      seed=seed)

res.glmnetclr$min.coef = Reduce(`+`,lapply(glmnetclr.res,function(x) x$glmnet.fit$beta[,x$index[1]]))[names(res.glmnetclr$min)]/Reduce(`+`,lapply(lapply(glmnetclr.res, function(x) x$glmnet.fit$beta[,x$index[1]]), function(x) x != 0))[names(res.glmnetclr$min)]
res.glmnetclr$`1se.coef` = Reduce(`+`,lapply(glmnetclr.res,function(x) x$glmnet.fit$beta[,x$index[2]]))[names(res.glmnetclr$`1se`)]/Reduce(`+`,lapply(lapply(glmnetclr.res, function(x) x$glmnet.fit$beta[,x$index[2]]), function(x) x != 0))[names(res.glmnetclr$`1se`)]

save(res.FLORAL,
     res.glmnetlog,
     res.glmnetrel,
     res.glmnetclr,
     file="allosct/res_new.Rdata")


### Run zeroSum with 100 cross validations
### Note that the fold splits are the same as used in FLORAL for each run

zeroSum.res <- list()

for (i in 1:mcv){
  
  zeroSum.res[[i]] <- zeroSum(x=as.matrix(cbind(covmat,log1p(df_taxa_abs))),
                              y=as.matrix(Surv(df_meta$OStime_landmark, df_meta$OSevent)),
                              family="cox",nFold=5,foldid = FLORAL.res[[i]]$foldid, lambda = FLORAL.res[[1]]$lambda,
                              penalty.factor = c(rep(0,6),rep(1,105)),
                              zeroSum.weights =  c(rep(0,6),rep(1,105)),
                              intercept=FALSE)
  
}

res.zeroSum <- list(min=table(unlist(lapply(zeroSum.res,function(x) x$variables.names[which(as.matrix(x$coef[[x$lambdaMinIndex]]!=0))])))/mcv,
                    `1se`=table(unlist(lapply(zeroSum.res,function(x) x$variables.names[which(as.matrix(x$coef[[x$lambda1SEIndex]]!=0))])))/mcv,
                    mcv=mcv,
                    seed=seed)

res.zeroSum$min.coef = Reduce(`+`,lapply(zeroSum.res,function(x) as.vector(x$coef[[x$lambdaMinIndex]])[x$variables.names %in% names(res.zeroSum$min)]))/Reduce(`+`,lapply(zeroSum.res,function(x) as.vector(x$coef[[x$lambdaMinIndex]])[x$variables.names %in% names(res.zeroSum$min)] != 0))
res.zeroSum$`1se.coef` = Reduce(`+`,lapply(zeroSum.res,function(x) as.vector(x$coef[[x$lambda1SEIndex]])[x$variables.names %in% names(res.zeroSum$min)]))/Reduce(`+`,lapply(zeroSum.res,function(x) as.vector(x$coef[[x$lambda1SEIndex]])[x$variables.names %in% names(res.zeroSum$min)] != 0))

save(res.FLORAL,
     res.zeroSum,
     res.glmnetlog,
     res.glmnetrel,
     res.glmnetclr,
     file="peri_engraftment_cv.Rdata")