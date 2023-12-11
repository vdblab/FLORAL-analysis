options(echo=TRUE)
args = commandArgs(trailingOnly=T)

library(dplyr)
library(tidyverse)
library(lubridate)
library(phyloseq)
library(survival)
library(metagenomeSeq)
library(corncob)
library(phyloseq)
library(ANCOMBC)
library(LDM)
library(lefser)
library(readr)
library(compositions)
library(ALDEx2)
source("https://raw.githubusercontent.com/FrederickHuangLin/ANCOM/8abbfb67ed1cb9ae366ab6349e0ad8f017ca6512/scripts/ancom_v2.1.R")

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

# Create relative abundance and CLR-transformed taxa matrices
relabmat <- apply(df_taxa_abs,2,function(y) y/rowSums(df_taxa_abs))
clrmat <- matrix(clr(df_taxa_abs),nrow=nrow(df_taxa_abs),ncol=ncol(df_taxa_abs))
colnames(clrmat) <- colnames(df_taxa_abs)

### Initialize a data frame to save results
Selected <- data.frame()

### Create feature matrix x, and survival outcome object y
x=df_taxa_abs
y=Surv(df_meta$OStime_landmark, df_meta$OSevent)

### Filter out samples collected after the event
x=x[y[,1] > 0,]
y=y[y[,1] > 0,]

### LDM

cox.fit = coxph(y ~ 1)
resid.Martingale = residuals(cox.fit, type='martingale')
clin <- data.frame(outcome=resid.Martingale)
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

### LEfSe

se0 <- SummarizedExperiment(assays=SimpleList(counts=t(x)),
                            colData=data.frame(outcome=as.factor(y[,2])))

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

meta_data <- data.frame(outcome=y[,2],id=1:length(y))
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
  fit <- summary(lm(log(relabmat[,j]+1e-8)~y[,2]))
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
  wil <- wilcox.test(relabmat[y[,2]==1,j],relabmat[y[,2]==0,j])
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

conds <- y[,2]
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

### metagenomeSeq

tstart <- proc.time()
yn <- y
rownames(yn) <- rownames(x)
test_obj<- newMRexperiment(
  counts = t(x), 
  phenoData = AnnotatedDataFrame(data.frame(outcome=yn[,2]))
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
META = sample_data(data.frame(outcome=as.factor(y[,2])))
physeq = phyloseq(OTU, META)
# physeq = filter_taxa(physeq, function(x) sum(x > 0) > (0.05*length(x)), TRUE) #after filtering nothing changes

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
META = sample_data(data.frame(outcome=as.factor(y[,2])))
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

Selected$n <- length(y)
Selected$p <- ncol(x)

save(Selected,file="peri_engraftment_various_methods.Rdata")
