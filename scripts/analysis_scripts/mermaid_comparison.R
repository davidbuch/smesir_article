library(smesir)
library(data.table)
setwd("../MERMAID-master/")

#------------------------------------------

dt_state <- readRDS("data/US_processed/dt_states.rds")

all_states <- sort(unique(subset(dt_state, !is.na(region))$state))

states <- all_states[44] #[23]

# Load merged case count data
dt <- as.data.table(readRDS("data/US_processed/merged_final_US_S2S.rds"))

dt[,date := as.Date(date),]
data.table::setorder(dt, state,date)

# Create data set with key variables
odt <- dt[order(state,date)][,list(
  date, 
  incid = incid_S2S, 
  incid_w = incid_S2S, 
  pop_size = pop_size, 
  ntest_w = 1.00 + ntest_S2S
), by = list(state)]

odt <- odt[order(state,date)]

odt <- odt[state == states]

odt$case_counts <- round(odt$incid)

J <- nrow(odt)
N <- odt$pop_size[1]
T_1 <- 1#odt$date[1]
gamma <- 1/10 # 1/mean_recovery_time
psi <- rep(0.5,7)/14 # incidence event probabilities (daily data)

epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = 1/gamma, incidence_probabilities = psi, 
                   discount_period_length = 60)

TEST_LENGTH <- 50
rolling_endpoints <- seq(100, nrow(odt) - TEST_LENGTH, by = TEST_LENGTH)
### SET IT UP FOR MERMAID
### FIGURE COLORSCHEMES
### GET THE OTHER STUDIES RUNNING ON THE CLUSTER
### EXTEND ROLLING FIT COMPARISON TO PERIOD WITH VACCINATIONS
### GO WORK OUT
i <- 1
while(i <= length(rolling_endpoints)){
  print(paste0("Iteration ", i))
  sel_train <- 1:rolling_endpoints[i]
  sel_test <- (rolling_endpoints[i] + 1):(rolling_endpoints[i] + TEST_LENGTH)
  
  sdat <- data.frame(cases = odt$case_counts[sel_train])
  sfit <- smesir(cases ~ 1, data = sdat, epi_params = epi_params)
  
  worst_rhat <- max(sapply(sfit$mcmc_diagnostics, function(region_diagnostics) max(region_diagnostics[,"Rhat"])))
  if(worst_rhat > 1.1){
    print(paste0("MCMC Sampler Failed to Converge on iteration ", i, ". Repeating this step..."))
    next
  }
  
  sforecast <- smesir_forecast(TEST_LENGTH, sfit)
  event_confints <- t(apply(sforecast$event_samples, 1, quantile, probs = c(0.09, 0.45, 0.81)))
  
  sel_plotwindow <- (min(sel_test) - 50):max(sel_test)
  matplot(sel_train, event_confints[sel_train, ], 
          xlab = "Day", ylab = "Cases",
          xlim = c(1, 300), ylim = c(0, max(event_confints)),
          type="l", main=paste0("Rolling Fit ", i), col="navy", lty=1)
  matplot(sel_test, event_confints[sel_test, ],
          type="l", col="navy", lty=2, add=T)
  lines(sel_train, odt$case_counts[sel_train], lwd=2, lty=1)
  lines(sel_test, odt$case_counts[sel_test], lwd=2, lty=2)
  abline(v = sel_test[1], lty=3)
  
  i <- i + 1
}

train_period_start <- as.Date("2020-03-15")
train_period_end <- as.Date("2020-12-07")
test_period_start <- as.Date("2020-12-08")
test_period_end <- as.Date("2020-12-21")
train_date_range <- seq.Date(train_period_start, train_period_end, by="day")
test_date_range <- seq.Date(test_period_start, test_period_end, by="day")

dt_train <- odt[(odt$date >= train_period_start) & 
                  (odt$date <= train_period_end), ]
dt_test <- odt[(odt$date >= test_period_start) & 
                 (odt$date <= test_period_end), ]

sdat <- data.frame(cases = dt_train$case_counts)

begin_t <- Sys.time()
sfit <- smesir(cases ~ 1, data = sdat, epi_params = epi_params)
end_t <- Sys.time()
print(end_t - begin_t)

print(paste0("Worst Rhat: ", max(sapply(sfit$mcmc_diagnostics, function(region_diagnostics) max(region_diagnostics[,"Rhat"])))))

n_train <- length(train_date_range)
n_test <- length(test_date_range)

sforecast <- smesir_forecast(n_test, sfit)
event_confints <- t(apply(sforecast$event_samples, 1, quantile, probs = c(0.09, 0.45, 0.81)))

sel_plotwindow <- c(train_date_range[train_date_range > "2020-11-07"], test_date_range)
matplot(train_date_range, event_confints[1:n_train,], type ="l", xlab = "Date", ylab = "Cases",
        xlim = c(as.Date("2020-11-07"), max(test_date_range)), ylim = c(0, max(event_confints)), col="navy", lty=1)
lines(train_date_range, dt_train$incid)

matplot(test_date_range, event_confints[(n_train + 1):n_test, ],
        type="l", col="navy", lty=2, add=T)
lines(sel_train, odt$case_counts[sel_train], lwd=2, lty=1)
lines(test_date_range, odt$case_counts[sel_test], lwd=2, lty=2)
abline(v = sel_test[1], lty=3)

matplot(test_date_range, event_confints[(n_train+1):(n_train + n_test),], type ="l", xlab = "Date", ylab = "New Cases",
        xlim = c(as.Date("2020-11-07"), max(test_date_range)), ylim = c(0, 24000), col="navy", lty=1)
lines(c(train_date_range, test_date_range), c(dt_train$incid, dt_test$incid))
# lines(test_date_range, dt_test$incid, lwd=1, lty=1)
abline(v=test_date_range[1])


beta_samps <- exp(sfit$design_matrices[[1]] %*% t(sfit$samples$Xi)) / gamma
IIP_samps <- sfit$samples$IIP
DISP_samps <- sfit$samples$DISP

beta_cis <- apply(beta_samps, 1, quantile, probs=c(0.25,0.75,0.95))
plot(odt$date, beta_cis[1,],type="l", ylim=c(min(beta_cis),max(beta_cis)))
lines(odt$date, beta_cis[2,])
lines(odt$date, beta_cis[3,])
lines(odt$date, rowMeans(beta_samps),lty="dashed")

vaccinations <- rep(0, J)

get_events_sample <- function(s){
  expected_events <- solve_events(solve_infections(beta_samps[,s],
                                gamma, T_1,
                                IIP_samps[s], N, vaccinations), psi)
  pred <- rnbinom(length(expected_events), size = 1/DISP_samps[s], mu = expected_events)
  return(pred)
}
event_samps <- sapply(1:length(IIP_samps), get_events_sample)
event_cis <- apply(event_samps, 1, quantile, probs=c(0.025,0.975))

plot(odt$date, odt$incid, ylim = c(0,max(event_cis)), xlim = c(min(odt$date), max(odt$date + 10)),
     type="l", main="Observed Events and Inferred Events")
lines(odt$date, event_cis[1,], col = "red")
lines(odt$date, event_cis[2,])
lines(odt$date, event_cis[3,])
lines(sforecast$event_confints[,1], col = "red")

