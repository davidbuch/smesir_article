plot_confint <- function(bounds,density=NULL,angle=45,border=NULL,col=NA,lty = par("lty"),x=NULL){
  if(is.null(x)) x <- 1:nrow(bounds)
  polygon(c(x,rev(x)),c(bounds[,1],rev(bounds[,2])),density,angle,border,col,lty)
}

compute_min_vulnerable <- function(smesir_fit,k){
  if(ncol(smesir_fit$response_matrix) == 1){
    k <-  1 # it doesn't matter what they said, k must be 1
    xi_samps <- smesir_fit$samples$Xi
    iip_samps <- smesir_fit$samples$IIP
  }else{
    xi_samps <- smesir_fit$samples$Xi[,,k]
    iip_samps <- smesir_fit$samples$IIP[,k]
  }
  nsamps <- length(iip_samps)
  
  vaccinations <- smesir_fit$vaccinations
  
  Y <- matrix(smesir_fit$response_matrix[,k],ncol = 1)
  J <- nrow(Y)
  
  T_1 <- smesir_fit$epi_params$outbreak_times[k]
  N <- smesir_fit$epi_params$region_populations[k]
  psi <- smesir_fit$epi_params$incidence_probabilities
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  
  
  beta_samps <- smesir_fit$design_matrices[[k]] %*% t(xi_samps)
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  
  sus_samps <- matrix(nrow = J, ncol = nsamps)
  for(s in 1:nsamps){
    sus_samps[,s] <- solve_susceptible(beta_samps[,s],gamma,T_1,iip_samps[s],N, vaccinations[,k])
  }

  msus <- mean(sus_samps[J,])
  sesus <- 2*sd(sus_samps[J,])
  return(c(msus,msus-sesus, msus+sesus))
}


min_vulnerable <- function(smesir_fit,k){
  if(ncol(smesir_fit$response_matrix) == 1){
    k <-  1 # it doesn't matter what they said, k must be 1
    xi_samps <- smesir_fit$samples$Xi
    iip_samps <- smesir_fit$samples$IIP
  }else{
    xi_samps <- smesir_fit$samples$Xi[,,k]
    iip_samps <- smesir_fit$samples$IIP[,k]
  }
  nsamps <- length(iip_samps)
  
  vaccinations <- smesir_fit$vaccinations
  
  Y <- matrix(smesir_fit$response_matrix[,k],ncol = 1)
  J <- nrow(Y)
  
  T_1 <- smesir_fit$epi_params$outbreak_times[k]
  N <- smesir_fit$epi_params$region_populations[k]
  psi <- smesir_fit$epi_params$incidence_probabilities
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  
  
  beta_samps <- smesir_fit$design_matrices[[k]] %*% t(xi_samps)
  gamma <- 1/smesir_fit$epi_params$mean_removal_time

  sus_samps <- matrix(nrow = J, ncol = nsamps)
  for(s in 1:nsamps){
    sus_samps[,s] <- solve_susceptible(beta_samps[,s],gamma,T_1,iip_samps[s],N, vaccinations[,k])
  }
  print(paste0(smesir_fit$region_names[k],": ", round(mean(sus_samps[J,]),2), " +/- (",  round(2*sd(sus_samps[J,]),4), ")"))
  
  sus_expected <- rowMeans(sus_samps)
  sus_CI <- t(apply(sus_samps,1,function(x) quantile(x, c(0.025,0.975))))
  
  plot(sus_expected, type = "l", lty = "dashed",
       ylim = c(0, 1),
       ylab = "proportion vulnerable", xlab = "weeks",
       main = smesir_fit$region_names[k])
  plot_confint(sus_CI, density = 15, col = "blue")
}

plot_smoothing_dsn <- function(smesir_fit,k){
  if(ncol(smesir_fit$response_matrix) == 1){
    k <-  1 # it doesn't matter what they said, k must be 1
    xi_samps <- smesir_fit$samples$Xi
    iip_samps <- smesir_fit$samples$IIP
  }else{
    xi_samps <- smesir_fit$samples$Xi[,,k]
    iip_samps <- smesir_fit$samples$IIP[,k]

  }
  nsamps <- length(iip_samps)
  
  vaccinations <- smesir_fit$vaccinations
  
  Y <- matrix(smesir_fit$response_matrix[,k],ncol = 1)
  J <- nrow(Y)
  
  T_1 <- smesir_fit$epi_params$outbreak_times[k]
  N <- smesir_fit$epi_params$region_populations[k]
  psi <- smesir_fit$epi_params$incidence_probabilities
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  
  beta_samps <- array(dim = c(J,nsamps))
  deaths_samps <- array(dim = c(J,nsamps))
  deaths_CI <- array(dim = c(J,2))
  nu_samps <- array(dim = c(J,nsamps))
  nu_CI <- array(dim = c(J,2))
  e_deaths_comp <- function(beta_samp,iip_samp){
    solve_events(solve_infections(beta_samp, gamma, T_1, iip_samp, N, vaccinations[,k]),psi[,k])
  }
  infections_comp <- function(beta_samp,iip_samp){
    solve_infections(beta_samp, gamma, T_1, iip_samp, N, vaccinations[,k])
  }
  
  beta_samps <- smesir_fit$design_matrices[[k]] %*% t(xi_samps)
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  R0_expected <- rowMeans(beta_samps)/gamma
  R0_CI <- t(apply(beta_samps,1,function(x) quantile(x, c(0.025,0.975))))/gamma
  
  for(s in 1:nsamps){
    deaths_samps[,s] <- e_deaths_comp(beta_samps[,s],iip_samps[s])
    nu_samps[,s] <- solve_infections(beta_samps[,s],gamma,T_1,iip_samps[s],N,vaccinations[,k])
  }
  deaths_CI <- t(apply(deaths_samps,1,function(x) quantile(x, c(0.025,0.975))))
  nu_expected <- cumsum(rowMeans(nu_samps))/N
  nu_CI <- apply(t(apply(nu_samps,1,function(x) quantile(x, c(0.025,0.975)))),2,cumsum)/N
  
  plot(R0_expected, type = "l", lty = "dashed",
       ylim = c(0, max(R0_CI)),
       ylab = "transmission rate", xlab = "weeks",
       main = smesir_fit$region_names[k])
  plot_confint(R0_CI, density = 15, col = "blue")
  abline(h = 1, col = "red")

  plot(Y, ylab = "deaths", xlab = "weeks", type = "l",
       ylim = c(0, max(Y,deaths_CI)), main = smesir_fit$region_names[k])
  plot_confint(deaths_CI, density = 15, col = "blue")
}

plot_cumulative_infections <- function(smesir_fit,k){
  if(ncol(smesir_fit$response_matrix) == 1){
    k <-  1 # it doesn't matter what they said, k must be 1
    xi_samps <- smesir_fit$samples$Xi
    iip_samps <- smesir_fit$samples$IIP
  }else{
    xi_samps <- smesir_fit$samples$Xi[,,k]
    iip_samps <- smesir_fit$samples$IIP[s,k]
    
  }
  nsamps <- length(iip_samps)
  
  Y <- matrix(smesir_fit$response_matrix[,k],ncol = 1)
  J <- nrow(Y)
  
  T_1 <- smesir_fit$epi_params$outbreak_times[k]
  N <- smesir_fit$epi_params$region_populations[k]
  psi <- smesir_fit$epi_params$incidence_probabilities
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  
  beta_samps <- array(dim = c(J,nsamps))
  nu_samps <- array(dim = c(J,nsamps))
  nu_CI <- array(dim = c(J,2))
  infections_comp <- function(beta_samp,iip_samp){
    solve_infections(beta_samp, gamma, T_1, iip_samp, N, vaccinations[,k])
  }
  
  beta_samps <- smesir_fit$design_matrices[[k]] %*% t(xi_samps)
  gamma <- 1/smesir_fit$epi_params$mean_removal_time
  R0_expected <- rowMeans(beta_samps)/gamma
  R0_CI <- t(apply(beta_samps,1,function(x) quantile(x, c(0.025,0.975))))/gamma
  
  for(s in 1:nsamps){
    deaths_samps[,s] <- e_deaths_comp(beta_samps[,s],iip_samps)
    nu_samps[,s] <- solve_infections(beta_samp,gamma,T_1,iip_samp,N, vaccinations[,k])
  }
  deaths_CI <- t(apply(deaths_samps,1,function(x) quantile(x, c(0.025,0.975))))
  nu_expected <- cumsum(rowMeans(nu_samps))/N
  nu_CI <- t(apply(t(apply(nu_samps,1,function(x) quantile(x, c(0.025,0.975)))),2,cumsum))/N
  
  plot(R0_expected, type = "l", lty = "dashed", 
       ylim = c(0, max(R0_CI)),
       ylab = "transmission rate", xlab = "weeks",
       main = smesir_fit$region_names[k])
  plot_confint(R0_CI, density = 15, col = "blue")
  abline(h = 1, col = "red")
  
  plot(nu_expected, type = "l", lty = "dashed", 
       ylim = c(0, max(R0_CI)),
       ylab = "transmission rate", xlab = "weeks",
       main = smesir_fit$region_names[k])
  plot_confint(nu_CI, density = 15, col = "blue")
  abline(h = 1, col = "red")
  
  plot(Y, ylab = "deaths", xlab = "weeks", type = "l", 
       ylim = c(0, max(Y,deaths_CI)), main = smesir_fit$region_names[k])
  plot_confint(deaths_CI, density = 15, col = "blue")
}



plot_smoothing_dsn_old <- function(sfit){
  Y <- sfit$response_matrix
  J <- nrow(Y); K <- ncol(Y)
  
  T_1 <- sfit$epi_params$outbreak_times
  N <- sfit$epi_params$region_populations
  psi <- sfit$epi_params$incidence_probabilities
  gamma <- 1/sfit$epi_params$mean_removal_time
  
  if(K == 1){
    nsamps <- nrow(sfit$samples$Xi)
    
    beta_samps <- array(dim = c(J,nsamps))
    deaths_samps <- array(dim = c(J,nsamps))
    deaths_CI <- array(dim = c(J,2))
    
    e_deaths_comp <- function(beta_samp,iip_samp){
      solve_events(solve_infections(beta_samp, gamma, T_1, iip_samp, N),psi)
    }
    beta_samps <- sfit$design_matrices[[1]] %*% t(sfit$samples$Xi)
    for(s in 1:nsamps){
      deaths_samps[,s] <- e_deaths_comp(beta_samps[,s],sfit$samples$IIP[s])
    }
    deaths_CI <- t(apply(deaths_samps,1,function(x) quantile(x, c(0.025,0.975))))
    
    plot(Y, ylab = "deaths", xlab = "weeks", type = "l", 
         ylim = c(0, max(Y,deaths_CI)), main = sfit$region_names)
    plot_confint(deaths_CI, density = 15, col = "blue")
    
  }else{
    nsamps <- nrow(sfit$samples$Xi[,,1])
    
    beta_samps <- array(dim = c(J,nsamps,K))
    deaths_samps <- array(dim = c(J,nsamps,K))
    deaths_CI <- array(dim = c(J,2,K))
    for(k in 1:K){
      e_deaths_comp <- function(beta_samp,iip_samp){
        solve_events(solve_infections(beta_samp, gamma, T_1[k], iip_samp, N[k]),psi[,k])
      }
      beta_samps[,,k] <- sfit$design_matrices[[k]] %*% t(sfit$samples$Xi[,,k])
      for(s in 1:nsamps){
        deaths_samps[,s,k] <- e_deaths_comp(beta_samps[,s,k],sfit$samples$IIP[s,k])
      }
      deaths_CI[,,k] <- t(apply(deaths_samps[,,k],1,function(x) quantile(x, c(0.025,0.975))))
      
      plot(Y[,k], ylab = "deaths", xlab = "weeks", type = "l", 
           ylim = c(0, max(Y[,k],deaths_CI[,,k])), main = sfit$region_names[k])
      plot_confint(deaths_CI[,,k], density = 15, col = "blue")
    }
  }

}

