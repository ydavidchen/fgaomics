# Preliminary LASSO

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(caret)
library(glmnet)

MODE <- "classification"
SEP <- "=============="

load(paste0(DIR_OUT,"240217a_mega.RData"))
mega_expr <- winsorize(mega_expr, -5, 5)
mega_expr <- t(mega_expr)
identical(rownames(mega_expr), meta_clin$sample)

if(MODE == "classification") {
  meta_clin$Label <- 1 * meta_clin$highFga
  reg_fam <- "binomial"
  out_type <- "class"
} else {
  meta_clin$Label <- meta_clin$normFga
  reg_fam <- "gaussian"
  out_type <- "response"
}

wrapper <- function(X, y, train_size=0.8) {
  #'@description Wrapper to run experiment
  idx_train <- createDataPartition(y, p=train_size, list=FALSE)
  
  Xtrain <- X[idx_train, ]
  Xval <- X[-idx_train, ]
  
  ytrain <- y[idx_train]
  yval <- y[-idx_train]

  mod <- glmnet(Xtrain, ytrain, family=reg_fam, type.measure="mse")
  
  ypred_val <- predict(mod, Xval, s=min(mod$lambda), type=out_type)
  
  if(MODE == "classification") {
    print( confusionMatrix(factor(yval), factor(ypred_val)) )
  } else {
    print( eval_median_cut(ypred_val[,1], yval) )
    pROC::roc(yval>=median(yval), ypred_val, plot=TRUE)
  }
}

## Overall:
wrapper(mega_expr, meta_clin$Label)

## Individual macrogroups:
for(mg in MACRO) {
  print(paste(SEP, mg, SEP))
  
  clin_mg <- subset(meta_clin, Macro == mg)
  expr_mg <- mega_expr[rownames(mega_expr) %in% clin_mg$sample, ]
  
  stopifnot(identical(rownames(expr_mg), clin_mg$sample)) #required checkpoint
  
  wrapper(expr_mg, clin_mg$Label)
  
  print(paste0(SEP, SEP, SEP))
}
