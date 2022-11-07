library(smesir)
library(ggplot2)
library(latex2exp)
set.seed(123)

data("covid_cases")
data("covid_deaths")
data("covid_vaccinations")
data("delta_prevalence")
data("state_populations")
data("ifr_by_state")

us_states <- colnames(covid_cases)

# Set "outbreak times" equal to the first week of cases
T_1 <- apply(covid_cases, 2, function(cvec) match(TRUE, cvec > 0))

# Set "time to death" distribution based on clinical literature
time_to_death <- c(0.008, 0.145, 0.312, 0.283, 0.154, 0.064, 0.024, 0.006, 0.004)

# We try various combinations over the plausible ranges of the remaining fixed parameters
incidence_probabilities <- time_to_death%*%t(ifr_by_state[,2])
infectious_period <- 1.5
lengthscale <- 8

epi_params <- list(region_populations = state_populations, 
                   outbreak_times = T_1,
                   mean_removal_time = infectious_period, 
                   incidence_probabilities = incidence_probabilities, # median IFR guess
                  discount_period_length = 90,
                  discount_period_disp = 10)

state_prior <- list(V0 = c(10,10,0.000001),
                    IGSR = matrix(rep(c(2.01,0.101),3), nrow = 3, ncol = 2, byrow = TRUE),
                    expected_initial_infected_population = 10000,
                    expected_dispersion = 5*sqrt(2/3.14159), # pretty diffuse, N+ with sd = 5
                    ell = lengthscale)

outbreak_data <- list(deaths = covid_deaths)

state_fit <- smesir(deaths ~ 1, data = outbreak_data, 
                    epi_params = epi_params,
                    vaccinations = covid_vaccinations,
                    prior = state_prior, 
                    region_names = us_states)
min(sapply(state_fit$mcmc_diagnostics, function(di) min(di[,2])), na.rm = T)
max(sapply(state_fit$mcmc_diagnostics, function(di) max(di[,1])), na.rm = T)

betas <- betas_lower <- betas_upper <- matrix(0, nrow = nrow(delta_prevalence), ncol = ncol(delta_prevalence))
for(k in 1:length(us_states)){
  state_beta_samps <- state_fit$design_matrices[[k]] %*% t(state_fit$samples$Xi[,,k])
  betas[,k] <- rowMeans(state_beta_samps)
}


sel <- delta_prevalence > 0.0 & delta_prevalence < 1.0
ggplot(mapping = aes(x = delta_prevalence[sel], y=betas[sel])) + 
  geom_point() + 
  geom_smooth(method = "loess", formula = 'y ~ x') +
  labs(title = TeX("Nonparametric $\\widehat{\\log\\beta}_{t,k}$ vs. Delta Prevalence $X_{t,k}$"),
       x = "Delta Prevalence as Percent of New Cases",
       y = TeX("$\\widehat{\\log\\beta}_{t,k}$")) + 
  theme(text = element_text(size=15))
ggsave(filename = "beta_vs_delta.png", path = "output/us_analysis", 
       width = 7, height = 7, device='png', dpi=300)


