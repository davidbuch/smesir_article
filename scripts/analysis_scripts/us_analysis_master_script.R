library(smesir)

slurm_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

endpoint <- c(74, 77, 90)[ceiling(slurm_task_id/27)] # weeks of data
version <- ((slurm_task_id - 1) %% 27) + 1 # used later to select epi parameters

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
param_options <- expand.grid(ifr_level = 1:3, # index over using lower, middle, uppper end of range
                             infectious_period = c(1,1.5,2.0),
                             lengthscale = c(8,12,16)
)

ifr_level <- param_options$ifr_level[version]
infectious_period <- param_options$infectious_period[version]
lengthscale <- param_options$lengthscale[version]

epi_params <- list(region_populations = state_populations, 
                   outbreak_times = T_1,
                   mean_removal_time = infectious_period, 
                   incidence_probabilities = 
                     time_to_death%*%t(ifr_by_state[,ifr_level]),
                   discount_period_length = 90,
                   discount_period_disp = 10)

state_prior <- list(V0 = c(10,10,0.000001),
                    IGSR = matrix(rep(c(2.01,0.101),3), nrow = 3, ncol = 2, byrow = TRUE),
                    expected_initial_infected_population = 10000,
                    expected_dispersion = 5*sqrt(2/3.14159), # pretty diffuse, N+ with sd = 5
                    ell = lengthscale)

outbreak_data <- list(deaths = covid_deaths[1:endpoint,],
                      delta = delta_prevalence[1:endpoint,],
                      delta_sq = delta_prevalence[1:endpoint,]^2)

if(endpoint == 74){ # endpoint 74, no delta covariate
  state_fit <- smesir(deaths ~ 1, data = outbreak_data, 
                      epi_params = epi_params,
                      vaccinations = covid_vaccinations[1:endpoint,],
                      prior = state_prior, 
                      region_names = us_states)
}else{
  state_fit <- smesir(deaths ~ delta_sq, data = outbreak_data, 
                      epi_params = epi_params,
                      vaccinations = covid_vaccinations[1:endpoint,], 
                      prior = state_prior, 
                      region_names = us_states)
}

saveRDS(state_fit, file = paste0("./output/us_analysis/sfit_",endpoint,"_",version,".rds"))