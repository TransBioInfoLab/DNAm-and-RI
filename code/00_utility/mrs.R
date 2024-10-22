# Auxilary function
# when a = NULL, the function conducts elastic net
# when a = 0, the function conducts ridge regression
# when a = 1, the function conducts lasso
library(survival)
library(survminer)
library(glmnet)
library(doParallel)
library(foreach)
fit_glmnet <- function(cpgs, beta, y, k = 10, a = NULL, seed = 1080,
                       fixed_df = NULL, scale = F, center = F,
                       lambda = NULL,
                       type = "regression",...){
  
  x <- t(beta[cpgs,])
  
  if(!is.null(fixed_df)) {
    pf <- c(rep(1, ncol(x)), rep(0, ncol(fixed_df)))
    x <- cbind(x, fixed_df)
  } else {
    pf <- rep(1, ncol(x))
  }
  
  # remove NA as glmnet cannot handle missing value
  x <- x[!is.na(y),]
  y <- y[!is.na(y)]
  
  if(scale) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)])
  }
  if(center) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)], scale = F)
  }
  x <- as.matrix(x)
  #y <- scale(y)
  
  if(is.null(a)){
    a <- seq(0.1, 1, 0.1)
    # cross-validation
    doParallel::registerDoParallel(4)
    search <- foreach(i = a, .combine = rbind) %do% {
      set.seed(seed)
      foldid <- sample(1:k,size=length(y),replace=TRUE)
      cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                      foldid = foldid, alpha = i, penalty.factor = pf, ...)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    print(search)
    if(type == "survival"){
      cv <- search[which.max(search$cvm), ]
      cv_min <- cv$cvm
    } else {
      cv <- search[which.min(search$cvm), ]
      cv_min <- cv$cvm
    }
    l <- cv$lambda.1se
    alpha <- cv$alpha
  } else {
    set.seed(seed)
    foldid <- sample(1:k,size=length(y),replace=TRUE)
    cv <- cv.glmnet(x = x, y = y, lambda = lambda, 
                    foldid = foldid, alpha = a, penalty.factor = pf,...)
    l <- cv$lambda.1se
    alpha <- a

    cv_min <- cv$cvm[cv$lambda == l]

  }
  
  # Fit model
  set.seed(seed)
  ridge <- glmnet(x = x, 
                  y = y, 
                  lambda = l, 
                  alpha = alpha,
                  penalty.factor = pf,
                  ...)
  
  if(type == "survival"){
    coeff <- as.numeric(coef(ridge))
    names(coeff) <- colnames(x)
  } else {
    coeff <- as.numeric(coef(ridge))
    names(coeff) <- c("intersect", colnames(x))
  }
  
  
  mod <- list(
    coeff = coeff,
    cv_min = cv_min,
    lambda_min = l,
    a = alpha,
    mod = ridge
  )
  
  return(mod)
  
}
# MRS 
get_MRS <- function(cpg.weight, betaset, scale = F, center = F){
  
  ## Match CpGs 
  common.cpg <- intersect(
    names(cpg.weight), rownames(betaset)
  )
  
  beta.common <- betaset[common.cpg,]
  cpg.weight.common <- cpg.weight[common.cpg]
  
  # ## Scaling CpGs
  if(scale) {
    beta.common <- t(scale(t(beta.common)))
  }
  if(center) {
    beta.common <- t(scale(t(beta.common), scale = F))
  }
  
  MRS <- (cpg.weight.common %*% beta.common)[1,]
  
  return(
    list(
      MRS = MRS,
      common.cpg = common.cpg
    )
  )
}
# KM plot
KM_plot <- function(test_var, time_var, event_var, pheno_mat, cut = "median", conf.int = F, palette = "jco", ...){
  
  if(is.numeric(test_var)){
    
    if(cut == "median"){
      m <- median(test_var)
    }
    
    if(cut == "mean"){
      m <- mean(test_var)
    }
    
    if(cut == "maxstat"){
      
      df <- data.frame(cluster = test_var, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
      m <- maxstat::maxstat.test(Surv(time, death) ~ cluster, data = df, smethod = "LogRank")$estimate
      
    }
    gene_cut <- ifelse(test_var < m, "low", "high")
    
  } else {
    
    gene_cut = test_var
    
  }
  
  df <- data.frame(Resilience = gene_cut, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
  
  fo <- as.formula(paste0("Surv(time, death) ~ Resilience"))
  fit <- surv_fit(fo, data = df)
  p <- round(-log10(survdiff(fo, data = df)$p),4)
  
  survminer::ggsurvplot(
    fit,
    data = df,
    size = 1,
    conf.int = conf.int,
    pval = T,
    risk.table = TRUE,
    palette = palette,
    ...
  )
}
