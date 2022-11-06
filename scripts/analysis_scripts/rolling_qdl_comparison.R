library(smesir)
library(data.table)

#------------------------------------------

dt_state <- readRDS("data/data_from_qdl2021/dt_states.rds")

all_states <- sort(unique(subset(dt_state, !is.na(region))$state))

states <- all_states[44] #[23]

# Load merged case count data
dt <- as.data.table(readRDS("data/data_from_qdl2021/merged_final_US_S2S.rds"))

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
T_1 <- 1
gamma <- 1/10 # 1/mean_recovery_time
psi <- rep(0.5,7)/14 # incidence event probabilities (daily data)

epi_params <- list(region_populations = N, outbreak_times = T_1,
                   mean_removal_time = 1/gamma, incidence_probabilities = psi, 
                   discount_period_length = 60)

TEST_LENGTH <- 50
rolling_endpoints <- seq(100, nrow(odt) - TEST_LENGTH, by = TEST_LENGTH)
monthly <- seq.Date(as.Date("2020-03-01"), as.Date("2020-12-01"), by = "month")

png("output/us_analysis/cases/rolling_fit.png", width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
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
  event_confints <- t(apply(sforecast$event_samples, 1, quantile, probs = c(0.1, 0.5, 0.9)))
  
  sel_plotwindow <- (min(sel_test) - 50):max(sel_test)
  matplot(odt$date[sel_train], event_confints[sel_train, ], 
          xlab = "Day", ylab = "New Cases",
          xlim = c(odt$date[1], odt$date[300]), ylim = c(0,50000), #ylim = c(0, max(event_confints)),
          type="l", main=paste0("Rolling Fit ", i), col="navy", lty=1,
          xaxt='n')
  axis(1, monthly, format(monthly, "%b %y"), cex.axis = 1)
  matplot(odt$date[sel_test], event_confints[sel_test, ],
          type="l", col="navy", lty=2, add=T)
  lines(odt$date[sel_train], odt$case_counts[sel_train], lwd=2, lty=1)
  lines(odt$date[sel_test], odt$case_counts[sel_test], lwd=2, lty=2)
  abline(v = odt$date[sel_test[1]], lty=3)
  
  i <- i + 1
}
par(mfrow = c(1,1))
dev.off()

## QDL Comparison
train_period_start <- as.Date("2020-03-15")
train_period_end <- as.Date("2020-12-07")
test_period_start <- as.Date("2020-12-08")
test_period_end <- as.Date("2020-12-21")

train_date_range <- seq.Date(train_period_start, train_period_end, by="day")
test_date_range <- seq.Date(test_period_start, test_period_end, by="day")

sdat <- data.frame(cases = odt[odt$date %in% train_date_range]$case_counts)

sfit <- smesir(cases ~ 1, data = sdat, epi_params = epi_params)
print(paste0("Worst Rhat: ", max(sfit$mcmc_diagnostics[[1]][,1])))

n_train <- length(train_date_range)
n_test <- length(test_date_range)

sforecast <- smesir_forecast(n_test, sfit)
event_confints <- t(apply(sforecast$event_samples, 1, quantile, probs = c(0.20, 0.5, 0.80)))


png("output/us_analysis/cases/qdl_comparison.png", width = 8, height = 6, units = 'in', res = 300)
plot_range <- seq.Date(as.Date("2020-11-07"), as.Date("2020-12-21"), by="day")
plot(date_range, odt[odt$date %in% plot_range]$case_counts, type = "l",
     ylim = c(0, 24000), lty=1,
     xlab = "Date", ylab = "New Cases")
matplot(test_date_range, event_confints[(n_train + 1):(n_train + n_test), ],
        type="l", col="blue", lty=c(3,1,3), add=T)
abline(v = test_date_range[1], lty=1)
legend("bottomleft", 
       legend = c("Confirmed cases (reported)", "Confirmed cases (predicted)"), 
       col = c('black', 'blue'), lty = c(1,1))
dev.off()
