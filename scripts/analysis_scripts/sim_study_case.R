library(smesir)
library(xtable)
library(tidyverse)
library(gridExtra)
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

# Plot covariates
covariate_list <- list(X1, X2, X3, U1)
png("output/sim_study_case/simstudy_covariates.png", width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
for(i in 1:4){
  matplot(covariate_list[[i]], xlab = "time", ylab = NA, type = "l", 
          cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, lwd = 3, main = c("(a)","(b)","(c)","(d)")[i])
  #abline(v = Tobs, col = "red", lwd = 2, lty = "dashed")
  legend("topright",c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"), bg = "white", lwd = 3, lty = 1:5, col = 1:6)
}
par(mfrow = c(1,1))
dev.off()

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

## Plot transmission rate
png("output/sim_study_case/simstudy_beta.png", width = 8, height = 8, units = 'in', res = 300)
matplot(beta, type = "l", xlab = "time", ylab = NA, 
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, lwd = 3, main = "Transmission Rate")
legend("topright", c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"), bg = "white", lwd = 3, lty = 1:5, col = 1:6)
dev.off()


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

## Plot "observed" data
png("output/sim_study_case/simstudy_observed.png", width = 8, height = 8, units = 'in', res = 300)
matplot(Y, type = "l", xlab = "time", ylab = NA, 
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, lwd = 3, main = "Incidence Events")
legend("topright", c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5"), bg = "white", lwd = 3, lty = 1:5, col = 1:6)
dev.off()

## -- THIRD, FIT THE MODEL -- ##
Jo <- J # maybe later we will want to check forecasting performance, too
sdat <- list(deaths = Y[1:Jo,], X1 = X1[1:Jo,], X2 = X2[1:Jo,], X3 = X3[1:Jo,])


# epi params all specified 'correctly'
epi_params <- list(region_populations = N, outbreak_times = T_1, mean_removal_time = 1/gamma, incidence_probabilities = psi, discount_period_length = 0)

# priors should match generating distribution, 
prior <- list(ell = J/5, V0 = c(10,10,1e-16), 
              expected_initial_infected_population = 10.0,
              expected_dispersion = 0.5*sqrt(2/3.14159), # sd 0.5
              IGSR = matrix(rep(c(2.01,0.101),3), nrow = 3, ncol = 2, byrow = TRUE))

start_time <- Sys.time()
sfit <- smesir(deaths ~ X1 + X2 + X3, data = sdat, vaccinations = vaccinations, epi_params = epi_params, prior = prior, region_names = paste0("R", 1:K), quiet = TRUE)
end_time <- Sys.time()
print(end_time - start_time)

# Write a LaTeX table with each region's Max Rhat and Min ESS
fileConn <- file("output/sim_study_case/mcmc_diagnostics_table.txt")
mcmc_diagnostics_table <- xtable(t(sapply(sfit$mcmc_diagnostics, function(x) c(Rhat = max(x[,1]), ESS = min(x[,2])))), digits = c(0,3,0))
writeLines(print(mcmc_diagnostics_table), fileConn)
close(fileConn)

png("output/sim_study_case/simstudy_posteriors.png", width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(3,3))
prior_rgb <- col2rgb('moccasin') / 255
post_rgb <- col2rgb('dodgerblue2') / 255
truth_color <- 'salmon4'
# 9 plots (8 histograms, one time series with CI) global intercept, global effects 1, 2, 3, local intercept variance, local effect variance, dispersion parameter (region 1), iip (region 1), temporal random effect (region 1)
# 1. global intercept
global_intercept_prior <-  rnorm(10000, mean = 0, sd = sqrt(sfit$prior$V0[1]))
hist(global_intercept_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, main = "(a)", ylim = c(0,3.5), xlim = c(-1,1), xlab = NA)
abline(v = alpha0, col = truth_color, lwd = 3)
hist(sfit$samples$Xi0[,1], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, add = TRUE)

# 2. global effect 1
global_coef_prior <-  rnorm(10000, mean = 0, sd = sqrt(sfit$prior$V0[2]))
hist(global_coef_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, main = "(b)", ylim = c(0,3.5), xlim = c(-1,1), xlab = NA)
abline(v = phi0[1], col = truth_color, lwd = 3)
hist(sfit$samples$Xi0[,2], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, add = TRUE)

# 3. global effect 2
hist(global_coef_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, main = "(c)", ylim = c(0,3.5), xlim = c(-1,1), xlab = NA)
abline(v = phi0[2], col = truth_color, lwd = 3)
hist(sfit$samples$Xi0[,3], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, add = TRUE)

# 4. global effect 3
hist(global_coef_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, main = "(d)", ylim = c(0,3.5), xlim = c(-1,1), xlab = NA)
abline(v = phi0[3], col = truth_color, lwd = 3)
hist(sfit$samples$Xi0[,4], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(-20,20,0.1), freq = F, add = TRUE)

# 5. local intercept variance
local_int_var_prior <- 1/rgamma(10000, shape = sfit$prior$IGSR[1,1], rate = sfit$prior$IGSR[1,2])
hist(local_int_var_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = c(seq(0,50,0.005), Inf), freq = F, main = "(e)", ylim = c(0,35), xlim = c(0,0.1), xlab = NA)
abline(v = local_intercept_var, col = truth_color, lwd = 3)
hist(sfit$samples$V[,1], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(0,50,0.005), freq = F, add = TRUE)


# 6. local coef variance
local_coef_var_prior <- 1/rgamma(10000, shape = sfit$prior$IGSR[2,1], rate = sfit$prior$IGSR[2,2])
hist(local_coef_var_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = c(seq(0,50,0.005), Inf), freq = F, main = "(f)", ylim = c(0,35), xlim = c(0,0.1), xlab = NA)
abline(v = local_coef_var, col = truth_color, lwd = 3)
hist(sfit$samples$V[,2], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(0,50,0.005), freq = F, add = TRUE)


# 7. dispersion parameter (R1)
local_dispersion_prior <- abs(rnorm(10000, sd = sfit$prior$expected_dispersion/sqrt(2/3.14159)))
hist(local_dispersion_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(0,20,0.1), freq = F, main = "(g)", ylim = c(0,3), xlim = c(0,3), xlab = NA)
abline(v = dispersion[1], col = truth_color, lwd = 3)
hist(sfit$samples$DISP[,1], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(0,20,0.1), freq = F, add = TRUE)

# 8. initial infectious population (R1)
iip_prior <- rexp(10000, 1/10)
hist(iip_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(0,200,2), freq = F, main = "(h)", ylim = c(0,0.25), xlim = c(0,40), xlab = NA)
abline(v = IIP[1], col = truth_color, lwd = 3)
hist(sfit$samples$IIP[,1], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(0,200,2), freq = F, add = TRUE)

# 9. temporal random effect/unobserved heterogeneity (R3) 
re_posterior_confints <- t(apply(exp(sfit$samples$Xi[,5:11,3]%*%t(sfit$design_matrices$R3[,5:11])),2, function(x) quantile(x,c(0.025,0.975))))
plot(x = 1:J, y = exp(U1[,3]*phiu[3]), type = "l", col = truth_color, lwd = 3, xlab = NA, ylab = NA, main = "(i)", ylim = c(min(re_posterior_confints), max(re_posterior_confints)))
plot_confint <- function(bounds,density=NULL,angle=45,border=NULL,col=NA,lty = par("lty"),x=NULL){
  if(is.null(x)) x <- 1:nrow(bounds)
  polygon(c(x,rev(x)),c(bounds[,1],rev(bounds[,2])),density,angle,border,col,lty)
}
plot_confint(re_posterior_confints, density = 20, col = rgb(post_rgb[1], post_rgb[2], post_rgb[3]))
dev.off()


# re_posterior_confints <- t(apply(exp(sfit$samples$Xi[,,3]%*%t(sfit$design_matrices$R3)),2, function(x) quantile(x,c(0.025,0.975))))
# plot(x = 1:J, y = beta[,3], type = "l", col = truth_color, lwd = 2, xlab = NA, ylab = NA, main = "(i)", ylim = c(min(re_posterior_confints), max(re_posterior_confints)))
# plot_confint <- function(bounds,density=NULL,angle=45,border=NULL,col=NA,lty = par("lty"),x=NULL){
#   if(is.null(x)) x <- 1:nrow(bounds)
#   polygon(c(x,rev(x)),c(bounds[,1],rev(bounds[,2])),density,angle,border,col,lty)
# }
# plot_confint(re_posterior_confints, density = 20, col = rgb(post_rgb[1], post_rgb[2], post_rgb[3]))


## -- NOW, FIT ALL SINGLE REGION AND MULTI-REGION MODELS -- 
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


