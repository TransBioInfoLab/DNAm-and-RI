#######################################################################################################
# =================================================================================================== #
# Single CpG assocation test
# =================================================================================================== #
#######################################################################################################
lm_test <- function(beta, pheno, test_var, covariates = NULL, 
                    as_response = T, convert_to_M = T,
                    parallel = T,
                    robust = F,
                    scale = T,
                    fdr_method = "fdr", 
                    save = F, 
                    dir.save = NULL,
                    cores = 10,
                    prefix = "Framingham", ...) {
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  if(as_response) {
    res_var <- "cpg"
    test_var <- test_var
  } else {
    res_var <- test_var
    test_var <- "cpg"
  }
  
  if(convert_to_M) {
    mat <- minfi::logit2(beta)
  } else {
    mat <- beta
  }
  
  results <- plyr::adply(
    mat,
    .margins = 1,
    .fun = function(cg){
      if(scale) cg <- scale(cg)
      data <- data.frame(cpg = cg, pheno)
      suppressWarnings({
        get_lm_coef(
          res_var = res_var, 
          test_var = test_var,
          data = data,
          robust = robust,
          covariates = covariates,
          ...
        )
      })
      
    }, .id = "cpg", .parallel = parallel
  )
  
  colnames(results)[1] <- "cpg"
  
  if(robust) {
    degrees.freedom.value <- nrow(pheno[!is.na(pheno[[test_var]]),]) - (length(covariates) + 2)
    results$pval <- 2 * (1 - pt( abs(results$t_value), df = degrees.freedom.value) )
  }
  # Add fdr
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_model_results.csv"))
    )
  }
  return(results)
}
get_lm_coef <- function(res_var, test_var, data, covariates = NULL, robust = F, ...) {
  
  # make linear regression formula
  fo <-  paste0(res_var, " ~ ", test_var)
  if(!is.null(covariates)){
    fo <- paste0(fo, " + ", paste0(covariates, collapse = "+"))
  }
  formula <- stats::as.formula(fo)
  
  tryCatch({
    # fit regression
    if(robust) {
      lm_mod <- MASS::rlm(
        formula, 
        data = data,
        psi = MASS::psi.bisquare,
        maxit = 100
      )
    } else {
      lm_mod <- stats::lm(
        formula,
        data = data,
        ...
      )
    }
    
    # get summary and coefficient
    lm_coef <- data.frame(summary(lm_mod)$coefficients)
    
    coef_df <- lm_coef[grepl(test_var, rownames(lm_coef)),]
    
    # use clean name
    coef_df <- janitor::clean_names(coef_df)
    return(coef_df)
  }, error = function(e){
    # If error shows up
    return(NULL)
  })
  
}
# Add fdr
add_fdr <- function(results, method = "fdr"){
  
  pr <- results[[grep("^p",colnames(results), value = T)]]
  results$fdr <- p.adjust(pr, method = method)
  
  return(results)
  
}
# Bacon
bacon_adj <- function(data, 
                      est_var, 
                      z_var, 
                      std_var){
  
  ### 1. Compute genomic inflation factor before bacon adjustment
  data <- data %>% mutate(
    chisq = get(z_var)^2
  ) 
  
  # inflation factor - last term is median from chisq distrn with 1 df  
  inflationFactor <- median(data$chisq,na.rm = TRUE) / qchisq(0.5, 1)
  print("lambda")
  print(inflationFactor)
  
  ### 2. bacon analysis
  est <- data[[est_var]]
  se <- data[[std_var]]
  
  
  bc <- bacon(
    teststatistics = NULL,
    effectsizes = est,
    standarderrors = se,
    na.exclude = TRUE,
    verbose = F
  )
  #  posteriors(bc)
  # inflation factor
  print("lambda.bacon")
  print(inflation(bc))
  # bias
  print("estimate bias")
  print(bias(bc))
  print("estimates")
  print(estimates(bc))
  
  ### 3. Create final dataset
  data.with.inflation <- data.frame(
    data,
    Estimate.bacon = bacon::es(bc),
    StdErr.bacon = bacon::se(bc),
    pValue.bacon = pval(bc),
    fdr.bacon = p.adjust(pval(bc), method = "fdr"),
    stringsAsFactors = FALSE
  )
  
  print("o After bacon correction")
  print("Conventional lambda")
  lambda.con <- median((data.with.inflation$Estimate.bacon/data.with.inflation$StdErr.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
  print(lambda.con)
  
  # percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
  # percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
  # percent_2  <- 100 - percent_null - percent_1  
  bc2 <- bacon(
    teststatistics = NULL,
    effectsizes =  data.with.inflation$Estimate.bacon,
    standarderrors = data.with.inflation$StdErr.bacon,
    na.exclude = TRUE,
    priors = list(
      sigma = list(alpha = 1.28,  beta = 0.36), 
      mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
      epsilon = list(gamma = c(99, .5, .5)))
  )
  
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(estimates(bc2))
  data.with.inflation$chisq <- NULL
  
  inflation.stat <- data.frame(
    "Inflation.org" = inflationFactor,
    "Inflation.bacon" = inflation(bc),
    "Bias.bacon" = bias(bc),
    "Inflation.after.correction" = lambda.con,
    "Inflation.bacon.after.correction" = inflation(bc2),
    "Bias.bacon.after.correction" = bias(bc2)
  )
  
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "bacon.obj" = bc,
      "inflation.stat" = inflation.stat
    )
  )
}
# MRS 
get_MRS <- function(cpg.weight, betaset, scale = F){
  
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
  
  MRS <- (cpg.weight.common %*% beta.common)[1,]
  
  return(
    list(
      MRS = MRS,
      common.cpg = common.cpg
    )
  )
}
# ===================================================================================================
# Adjust M values function
# ===================================================================================================
methyl_adj <- function(mat, 
                       pheno, 
                       adjust_var, 
                       convert_to_M = T, 
                       return_to_beta = T,
                       parallel = T){
  
  if(convert_to_M){
    # transform beta to M
    mat <- minfi::logit2(mat)
  }
  mat <- as.matrix(mat)
  if(parallel) doParallel::registerDoParallel(10)
  
  resid_mat <- plyr::aaply(
    mat,
    .margins = 1,
    .fun = function(m){
      
      dat <- data.frame(M = m, pheno)
      # Create formula
      fo <- paste0("M ~ ", paste0(adjust_var, collapse = "+"))
      
      # Fit LM model
      lm_mod <- lm(as.formula(fo), data = dat)
      
      # Get residuals
      r <- resid(lm_mod)
      
      return(r)
    },.parallel = parallel
  )
  
  if(return_to_beta){
    # Transform back to beta
    resid_mat <- minfi::ilogit2(resid_mat)
  }
  
  return(resid_mat)
  
}
