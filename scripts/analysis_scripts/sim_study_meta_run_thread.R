# Create Containers
local_param_coverage <- array(dim = c(6,K,Replicates))
local_param_true_values <- array(dim = c(6,K,Replicates))
global_param_coverage <- array(dim = c(8,Replicates))
global_param_true_values <- array(dim = c(8,Replicates))

min_ess_vals <- array(dim = c(Replicates,6))
max_rhat <- array(dim = c(Replicates,6))

fit_summaries <- list()

# For each iteration, simulate parameters, data, and 
for(rpl in 1:Replicates){
  #print(paste("Replicate",rpl,"of",Replicates))
  
  ## Randomize true parameters 
  # Sample from the prior distributions
  local_intercept_var <- 1/rgamma(1,prior$IGSR[1,1],prior$IGSR[1,2])
  local_coef_var <- 1/rgamma(1,prior$IGSR[2,1],prior$IGSR[2,2])
  IIP <- 1 + rexp(K, 1/prior$expected_initial_infected_population)
  dispersion <- abs(rnorm(K, sd = 0.5))
  
  # simulate random global and local coeffs - (do so until we get valid beta vectors )
  beta_is_okay <- FALSE
  while(!beta_is_okay){
    alpha0 <- rnorm(1, mean = 2*gamma, sd = 0.4*gamma) # baseline R0 will be between 1.2 and 2.8
    phi0 <-  rnorm(3, mean = c(0,-0.5,0.5), sd = 0.2)
    alpha <- rnorm(K, alpha0, sd = sqrt(local_intercept_var))
    phi <- matrix(rnorm(length(phi0)*K, phi0, sd = sqrt(local_coef_var)), nrow = length(phi0), ncol = K)
    phiu <- rnorm(K, sd = 0.3)
    ## Compute regional beta from covariates, covariate effects, and "random effects"
    beta <- matrix(NA,J,K)
    for(k in 1:K){
      beta[,k] <- alpha[k] + X1[,k]*phi[1,k] + X2[,k]*phi[2,k] + X3[,k]*phi[3,k] + U1[,k]*phiu[k]
    }
    
    ## Simulate the incidence counts
    Y <- matrix(nrow = J, ncol = K)
    EY <- matrix(nrow = J, ncol = K)
    for(k in 1:K){
      EY[,k] <- solve_events(solve_infections(beta[,k],
                                    gamma, T_1[k],
                                    IIP[k], N[k],vaccinations),psi)
      Y[,k] <- rnbinom(J,size = 1/dispersion[k], mu = EY[,k])
    }
    # matplot(cbind(Y,EY),type = "l", main = dispersion)
    # beta must be nonnegative and we require regions to have at least a few days with events
    if(all(beta > 0) && apply(Y > 0,2,function(col) sum(col) >= 3)) beta_is_okay <- TRUE
  }

  sdat <- list(deaths = Y, X1 = X1, X2 = X2, X3 = X3)
  # Fit the model
  sfit <- smesir(deaths ~ X1 + X2 + X3, data = sdat, epi_params = epi_params,
                 region_names = paste0("R", 1:K), prior = prior)
  
  min_ess_vals[rpl,] <- sapply(sfit$mcmc_diagnostics,function(stats) min(stats[,2]))
  max_rhat[rpl,] <- sapply(sfit$mcmc_diagnostics,function(stats) max(stats[,1]))
  
  for(k in 1:K){
    localvals <- c(alpha[k],phi[,k],IIP[k],dispersion[k])
    local_param_coverage[,k,rpl] <- sfit$summary[[k]][,3] < localvals & sfit$summary[[k]][,4] > localvals
    local_param_true_values[,k,rpl] <- localvals
  }
  globalvals <- c(alpha0,phi0,local_intercept_var,local_coef_var,var(alpha - alpha0),var(c(phi-phi0)))
  global_param_coverage[,rpl] <- c(sfit$summary[[K+1]][1:6,3],sfit$summary[[K+1]][5:6,3]) < globalvals & c(sfit$summary[[K+1]][1:6,4],sfit$summary[[K+1]][5:6,4]) > globalvals
  global_param_true_values[,rpl] <- globalvals
  fit_summaries[[rpl]] <- sfit$summary
}

# local_coverage_rates <- apply(local_param_coverage,1,mean)
# names(local_coverage_rates) <- c("Int.", "X1", "X2", "X3", "IIP", "DISP")
# print(local_coverage_rates)
# 
# global_coverage_rates <- rowMeans(global_param_coverage)
# names(global_coverage_rates) <- c("Int.", "X1", "X2", "X3", "VLI", "VLC","VLI(S)", "VLC(S)")
# print(global_coverage_rates)

save(fit_summaries, min_ess_vals, max_rhat, local_param_coverage, global_param_coverage,
		   local_param_true_values, global_param_true_values,
		   file = paste0("output/sim_study_meta/", version,"_", thread,".RData"))