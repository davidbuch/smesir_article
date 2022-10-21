library(smesir)
data("covid_cases")
data("covid_deaths")
data("covid_vaccinations")
data("delta_prevalence")
data("state_populations")
data("ifr_by_state")
library(DescTools)

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

outbreak_data <- list(deaths = covid_deaths,
                      delta = delta_prevalence,
                      delta_sq = delta_prevalence^2)

state_fit <- smesir(deaths ~ delta_sq, data = outbreak_data, 
                    epi_params = epi_params,
                    vaccinations = covid_vaccinations,
                    prior = state_prior, 
                    region_names = us_states)
sapply(state_fit$mcmc_diagnostics, function(di) min(di[,2]))

local_effect_quantiles <- apply(state_fit$samples$Xi[,2,], 2, quantile, probs=c(0.025, 0.5, 0.975))
local_effect_quantiles <- cbind(local_effect_quantiles, global = quantile(state_fit$samples$Xi0[,2], probs = c(0.025, 0.5, 0.975)))

png("output/us_analysis/endpoints/usa_local_effects.png", width = 6, height = 12, units = 'in', res = 300)
PlotDot(x = local_effect_quantiles[2,], 
        args.errbars = list(from = local_effect_quantiles[1,], to = local_effect_quantiles[3,]), ,
        xlim = c(0,2),
        main = "Local Delta Variant Effects")
dev.off()


