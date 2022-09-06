library(smesir)
set.seed(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))

version_list <- c("true_epi", "misspecified_epi_1", "misspecified_epi_2", "radically_misspecified_epi")

## replicates per thread
Replicates <- 100

# Set Global Variables - Data dimensions and epidemiologic parameters
J <- 50; K <- 5
N <- c(2e5,6e5,5e5,3e5,2e5) # region populations
T_1 <- c(1,3,2,2,1) # outbreak times
gamma <- 2/3 # 1/mean_recovery_time
psi <- 0.01*c(0.25,0.5,0.25) # incidence event probabilities
vaccinations <- matrix(0, nrow = J, ncol = K)

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


## Begin Simulation Studies
thread <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# Use a common prior across all versions
prior <- list(ell = J/5, V0 = c(10,10,0.000001), 
              expected_initial_infected_population = 10.0,
              expected_dispersion = 0.5*sqrt(2/3.14159),
              IGSR = matrix(rep(c(2.01,0.101),3), nrow = 3, ncol = 2, byrow = TRUE))


version <- version_list[1] # "true_epi"
# epi_params will be used when fitting the model, but for simulation we use the original varables
epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = 1/gamma, incidence_probabilities = psi, discount_period_length = 0)
source("scripts/analysis_scripts/sim_study_meta_run_thread.R")

version <- version_list[2] # "misspecified_epi_1"
# epi_params will be used when fitting the model, but for simulation we use the original varables
epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = 1/gamma, incidence_probabilities = 0.9*psi, discount_period_length = 0)
source("scripts/analysis_scripts/sim_study_meta_run_thread.R")

version <- version_list[3] # "misspecified_epi_2"
# epi_params will be used when fitting the model, but for simulation we use the original varables
epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = (4/3)*1/gamma, incidence_probabilities = psi, discount_period_length = 0)
source("scripts/analysis_scripts/sim_study_meta_run_thread.R")

version <- version_list[4] # "radically_misspecified_epi"
# epi_params will be used when fitting the model, but for simulation we use the original varables
epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = (0.75)*1/gamma, incidence_probabilities = (1.5)*psi, discount_period_length = 0)
source("scripts/analysis_scripts/sim_study_meta_run_thread.R")
