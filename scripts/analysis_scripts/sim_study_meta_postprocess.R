library(xtable)
version_list <- c("true_epi", "misspecified_epi_1", "misspecified_epi_2", "radically_misspecified_epi")

FIT_SUMMARIES <- list()
MIN_ESS <- MAX_RHAT <- GLOBAL_COVERAGE <- GLOBAL <- NULL
LOCAL_COVERAGE <- LOCAL <- array(dim = c(6,5,1000))
for(version in version_list){
  for(thread in 1:10){
    load(paste0("output/sim_study_meta/", version,"_", thread,".RData"))
    
    if(is.null(MIN_ESS)){
      MIN_ESS <- min_ess_vals
      MAX_RHAT <- max_rhat
      GLOBAL_COVERAGE <- global_param_coverage
      GLOBAL <- global_param_true_values
      LOCAL_COVERAGE[,,1:1000] <- local_param_coverage
      LOCAL[,,1:1000] <- local_param_true_values
    }else{
      MIN_ESS <- rbind(MIN_ESS,min_ess_vals)
      MAX_RHAT <- rbind(MAX_RHAT,max_rhat)
      GLOBAL_COVERAGE <- cbind(GLOBAL_COVERAGE, global_param_coverage)
      GLOBAL <- cbind(GLOBAL, global_param_true_values)
      LOCAL_COVERAGE[,,1:100 + (thread-1)*100] <- local_param_coverage
      LOCAL[,,1:100 + (thread-1)*100] <- local_param_true_values
    }
    FIT_SUMMARIES <- append(FIT_SUMMARIES,fit_summaries)
  }
  
  # Make table about the posterior concentration
  # prior variance for the intercept and coefficients is 10
  # and the distribution is centered at zero
  global_posterior_sds <- matrix(nrow = 1000, ncol = 7)
  for(i in 1:length(FIT_SUMMARIES)){
    for(p in 1:7){ # There are 7 "global" params
      global_posterior_sds[i,p] <- FIT_SUMMARIES[[i]]$Global[p,2]
    }
  }
  global_posterior_sds <- colMeans(global_posterior_sds)
  prior_sds <- c(rep(sqrt(10),4),rep(.101^2/(((2.01 - 1)^2)*(2.01-2)),2),.2^2/(((3 - 1)^2)*(3 - 2)))
  
  conc_table <- rbind(prior_sds,global_posterior_sds)
  colnames(conc_table) <- names(FIT_SUMMARIES[[1]]$Global[1:7,2])
  rownames(conc_table) <- c("prior std. dev.", "posterior std. dev.")
  
  fileConn <- file(paste0("output/sim_study_meta/concentration_",version,".txt"))
  concentration_table <- xtable(t(conc_table))
  writeLines(print(concentration_table), fileConn)
  close(fileConn)
  
  # Make a table summarizing the coverage rate
  coverage_table <- matrix(nrow = 8, ncol = 6)
  colnames(coverage_table) <- c(paste0("R",1:5), "G")
  rownames(coverage_table) <- c("Intercept", "X1", "X2", "X3", "IIP","DISP","Local Int. Var.", "Local Coef. Var")
  coverage_table[1:6,1:5] <- apply(LOCAL_COVERAGE,c(1,2),mean) # local average coverage
  coverage_table[c(1:4,7:8),6] <- rowMeans(GLOBAL_COVERAGE[1:6,]) # global average coverage
  
  fileConn <- file(paste0("output/sim_study_meta/coverage_",version,".txt"))
  coverage_table <- xtable(coverage_table*100, digits = 0)
  writeLines(print(coverage_table), fileConn)
  close(fileConn)
}
