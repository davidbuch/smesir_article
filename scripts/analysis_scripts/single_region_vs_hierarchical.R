library(tidyverse)
library(gridExtra)
library(smesir)
set.seed(123)

J <- 50; K <- 5

## -- FIRST, GENERATE TRANSMISSION RATES -- ##

# Construct Covariates
# linear time (centered)
X1 <- scale(matrix(rep(1:J - J/2,K), nrow = J, ncol = K))

# lockdown
lockdown_start <- c(10,12,14,16,18)
lockdown_end <- c(20,22,24,26,28)
lockdown_restart <- c(40,38,36,34,32)
X2 <- matrix(0, nrow = J, ncol = K)
for(k in 1:K){
  X2[lockdown_start[k]:lockdown_end[k],k] <- 1
  X2[lockdown_restart[k]:J,k] <- 1
}

# holiday periods
X3 <- matrix(0, nrow = J, ncol = K)
X3[c(25,43:45),] <- 1

# sinusoidal covariate (unobserved)
period <- J/2 # period will be J/2
wn <- (2*3.14)/period
off <- (0:(K - 1))*period/K
U1 <- matrix(NA, nrow = J, ncol = K)
for(k in 1:K){
  U1[,k] <- sin(wn*(off[k] + 1:J)) # sin wave with period J/2
}


## Set Covariate Effects
alpha0 <- log(1.25)
phi0 <- log(c(1.00,0.75,1.25))

# scale of the local variation around the global effects
local_intercept_var <- 0.15^2
local_coef_var <- 0.15^2
alpha <- rnorm(K, alpha0, sd = sqrt(local_intercept_var))
phi <- matrix(rnorm(length(phi0)*K, phi0, sd = sqrt(local_intercept_var)), nrow = length(phi0), ncol = K)
phiu <- rnorm(K, sd = 0.3)

## Compute regional beta from covariates, covariate effects, and "random effects"
beta <- matrix(NA,J,K)
for(k in 1:K){
  beta[,k] <- exp(alpha[k] + X1[,k]*phi[1,k] + X2[,k]*phi[2,k] + X3[,k]*phi[3,k] + U1[,k]*phiu[k])
}

## -- SECOND, SIMULATE EVENT DATA -- ##

# set remaining parameters which will be "unknown"
dispersion <- abs(rnorm(K, sd = 0.5))
IIP <- 1 + rexp(K, 1/10)

# parameters which will be "known"
N <- c(2e5,6e5,5e5,3e5,2e5) # region populations
T_1 <- c(1,3,2,2,1) # outbreak times
gamma <- 2/3 # 1/mean_recovery_time
psi <- 0.01*c(0.25,0.5,0.25) # incidence event probabilities
vaccinations <- matrix(0, nrow = J, ncol = K)

Y <- matrix(nrow = J, ncol = K)
for(k in 1:K){
  Y[,k] <- rnbinom(J,size = 1/dispersion[k], mu = solve_events(solve_infections(beta[,k],gamma, T_1[k],IIP[k], N[k],vaccinations),psi))
}


# X <- cbind(12:J, X1[12:J,2], X2[12:J,2], X3[12:J,2])
# result_list <- BaySIR_MCMC(B = Y[12:J,2], I_D_0 = 18, 
#                            N = N[2], 
#                            X = X, 
#                            nu_alpha_1 = 300, 
#                            nu_alpha_2 = 200)
# result_list_2 <- BaySIR_MCMC(B = Y[12:J,2], I_D_0 = 18, 
#                            N = N[2], 
#                            X = X,
#                            nu_alpha_1 = 300, 
#                            nu_alpha_2 = 200)

prior <- list(ell = 15, V0 = c(10,10,1e-16), 
              expected_initial_infected_population = 10.0,
              expected_dispersion = 0.5*sqrt(2/3.14159), # sd 0.5
              IGSR = matrix(rep(c(2.01,0.101),3), nrow = 3, ncol = 2, byrow = TRUE))
epi_params <- list(region_populations = N, outbreak_times = T_1, mean_removal_time = 1/gamma, incidence_probabilities = psi, discount_period_length = 0)
sdat <- list(deaths = Y, X1 = X1, X2 = X2, X3 = X3)
sresG <- smesir(deaths ~ X1 + X2 + X3, data = sdat, epi_params = epi_params, prior = prior)

sres <- list()
for(k in 1:K){
  prior <- list(ell = 15, V0 = c(10,10), 
                expected_initial_infected_population = 10.0,
                expected_dispersion = 0.5*sqrt(2/3.14159), # sd 0.5
                IGSR = c(2.01,0.101))
  epi_params <- list(region_populations = N[k], outbreak_times = T_1[k], mean_removal_time = 1/gamma, incidence_probabilities = psi, discount_period_length = 0, lengthscale = 10)
  sdat <- data.frame(deaths = Y[,k], X1 = X1[,k], X2 = X2[,k], X3 = X3[,k])
  sres[[k]] <- smesir(deaths ~ X1 + X2 + X3, data = sdat, epi_params = epi_params, prior = prior)
}

intervalL <- list()
intervalG <- list()
for(k in 1:K){
  intervalL[[k]] <- cbind(
    t(exp(apply(sres[[k]]$design_matrices[[1]] %*% t(sres[[k]]$samples$Xi), 1, quantile, prob=c(0.025, 0.975)))),
    beta[,k]
  )
  intervalG[[k]] <- cbind(
    t(exp(apply(sresG$design_matrices[[k]] %*% t(sresG$samples$Xi[,,k]), 1, quantile, prob=c(0.025, 0.975)))),
    beta[,k]
  )
}
make_ribbonplot <- function(interval){
  ggplot() + 
    geom_ribbon(aes(x = 1:J, ymin = interval[,1], ymax = pmin(interval[,2],4))) + 
    ylim(0,4) + 
    labs(x = NULL, y = NULL) + 
    geom_line(aes(x = 1:J, y = interval[,3]), color = "white")
}

rplotsL <- lapply(intervalL, make_ribbonplot)
rplotsG <- lapply(intervalG, make_ribbonplot)
g1 <- grid.arrange(grobs = rplotsL,
             widths = c(1),
             layout_matrix = matrix(1:K, ncol = 1), 
             top = "Single-Region Model")
g2 <- grid.arrange(grobs = rplotsG,
                   widths = c(1),
                   layout_matrix = matrix(1:K, ncol = 1), 
                   top = "Hierarchical Model")
gB <- grid.arrange(grobs = list(g1,g2),
                   widths = c(1,1),
                   layout_matrix = matrix(1:2, nrow = 1), 
                   left = "Region",
                   top = "Dynamic Transmission Rates")

ydat <- data.frame(cbind(1:J,Y))
names(ydat) <- c("time", paste0('Y', 1:K))
ydat <- ydat %>% pivot_longer(paste0('Y', 1:K), names_to = "region", names_prefix = "Y", values_to = "count")
gT <- ggplot(ydat, aes(x = time, y = count, linetype = region)) + geom_line() + ggtitle("Regional Event Counts")

png("output/sim_study_case/sr_vs_hier_betas.png", width = 8, height = 8, units = 'in', res = 300)
g <- grid.arrange(grobs = list(gT, gB),
                  heights = c(2,5),
                  layout_matrix  = matrix(1:2, nrow = 2))
dev.off()

