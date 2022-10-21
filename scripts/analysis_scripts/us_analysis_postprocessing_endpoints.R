library(smesir)
library(lubridate)
library(DescTools)

# load covid data from the smesir package
data("covid_vaccinations")
data("covid_deaths")
data("delta_prevalence")
data("state_populations")
us_states <- names(state_populations)

# specify what models we are going to work with
version <- 14 # "default" epi parameters
endpoint <- c(74, 77, 90)
final_endpoint <- endpoint[length(endpoint)]
names(endpoint) <- as.character(endpoint)

# load previously fitted models and compute forecasts
sfit <- lapply(endpoint, function(e){
  readRDS(paste0("output/us_analysis/sfit_",e, "_", version, ".rds"))
})

sfor <- lapply(endpoint[1:2], function(e){
  smesir_forecast(final_endpoint - e, sfit[[as.character(e)]],
                  new_x = list(delta_sq = delta_prevalence[(e + 1):final_endpoint,]^2), 
                  new_vaccinations = covid_vaccinations[(e + 1):final_endpoint,],
                  fixed_dispersion = 10)
})

# prepare for plotting - load helper functions
source("scripts/analysis_scripts/plotting_functions.R")

# prepare for plotting - create monthly date vector
weeks <- ymd(rownames(covid_deaths))
yii <- year(min(weeks))
mii <- month(min(weeks))
monthly <- c()
while(yii < year(max(weeks)) || mii <= month(max(weeks))){
  monthly <- c(monthly, paste0(yii,"-",mii,"-01"))
  
  mii <- mii + 1
  if(mii > 12){
    yii <- yii + 1
    mii <- 1
  }
}
monthly <- ymd(monthly)


# Compare vulnerable population size at the end of 74 and 78 weeks

vulnerable74<- t(sfor[[names(endpoint[1])]]$vulnerable_confints[endpoint[1],,])
dimnames(vulnerable74) <- list(us_states,c("2.5%","50%","97.5%"))
order74 <- order(vulnerable74[,2])
vulnerable74 <- vulnerable74[order74,]

vulnerable78 <- t(sfor[[names(endpoint[1])]]$vulnerable_confints[78,,])
dimnames(vulnerable78) <- list(us_states,c("2.5%","50%","97.5%"))
vulnerable78 <- vulnerable78[order74,]

pdf("output/us_analysis/endpoints/remaining_vulnerable_endpoints.pdf", height = 11, width = 7)
PlotDot(x = vulnerable74[,2], args.errbars = list(from = vulnerable74[,1], to = vulnerable74[,3]))
PlotDot(x = vulnerable78[,2], args.errbars = list(from = vulnerable78[,1], to = vulnerable78[,3]), add = T)
dev.off()
# plot/compare forecasts for one state

## --- FORECAST COMPARISON ----
PLOTSTART <- 50 # lower bound on time

# Compare Forecasts U.S.A.
event_quantile1 <- rowSums(apply(sfor[[names(endpoint[1])]]$expected_samples,c(1,3),function(x) quantile(x, c(0.6))))
event_quantile3 <- rowSums(apply(sfor[[names(endpoint[2])]]$expected_samples,c(1,3),function(x) quantile(x, c(0.6))))

png("output/us_analysis/endpoints/usa_endpoint_forecasts.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:endpoint[3]], rowSums(covid_deaths[PLOTSTART:endpoint[3],]), xaxt = 'n', xlab = NA, ylab = NA, main = "US Weekly COVID-19 Deaths", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)

abline(v = weeks[endpoint[1]], col="moccasin")
abline(v = weeks[endpoint[2]], col ="dodgerblue2")
plot_confint(cbind(rowSums(sfor[[names(endpoint[1])]]$expected_confints[(endpoint[1] + 1):endpoint[3],1,]),rowSums(sfor[[names(endpoint[1])]]$expected_confints[(endpoint[1] + 1):endpoint[3],3,])),x=weeks[(endpoint[1] + 1):endpoint[3]],density=20,angle=90, col="moccasin")
plot_confint(cbind(rowSums(sfor[[names(endpoint[2])]]$expected_confints[(endpoint[2] + 1):endpoint[3],1,]),rowSums(sfor[[names(endpoint[2])]]$expected_confints[(endpoint[2] + 1):endpoint[3],3,])),x=weeks[(endpoint[2] + 1):endpoint[3]],density=20,angle=0, col="dodgerblue2")
lines(x=weeks[(endpoint[1] + 1):endpoint[3]],event_quantile1[(endpoint[1] + 1):endpoint[3]],lty = 3,lwd = 4, col="moccasin")
lines(x=weeks[(endpoint[2] + 1):endpoint[3]],event_quantile3[(endpoint[2] + 1):endpoint[3]],lty = 4,lwd = 4, col="dodgerblue2")
legend("topleft",legend = c("Forecast 95% CI after 74 weeks","Forecast 95% CI after 77 weeks","Forecast 0.60 Quantile after 74 weeks","Forecast 0.60 Quantile after 77 weeks"), fill = c("moccasin", "dodgerblue2", NA, NA), density = c(20,20,0,0), angle = c(90,0,NA,NA), border = c("moccasin", "dodgerblue2", NA, NA), lty = c(NA,NA,3,4), lwd = c(NA, NA, 4,4), col = c(NA,NA, "moccasin", "dodgerblue2"),cex = 1, bty = "o", bg = "white" )
dev.off()

# Compare Forecasts California
k <- which(us_states == "california")
cal_event_quantile1 <- apply(sfor[[names(endpoint[1])]]$event_samples[,,k],1,function(x) quantile(x, c(0.85)))
cal_event_quantile3 <- apply(sfor[[names(endpoint[2])]]$event_samples[,,k],1,function(x) quantile(x, c(0.85)))

png("output/us_analysis/endpoints/cal_endpoint_forecasts.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:endpoint[3]], covid_deaths[PLOTSTART:endpoint[3],k], xaxt = 'n', xlab = NA, ylab = NA, main = "CA Weekly COVID-19 Deaths", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)
abline(v = weeks[endpoint[1]], col="moccasin")
abline(v = weeks[endpoint[2]], col = "dodgerblue2")
plot_confint(sfor[[names(endpoint[1])]]$expected_confints[(endpoint[1] + 1):endpoint[3],c(1,3),k],x=weeks[(endpoint[1] + 1):endpoint[3]],density=30,angle=90,col="moccasin")
plot_confint(sfor[[names(endpoint[2])]]$expected_confints[(endpoint[2] + 1):endpoint[3],c(1,3),k],x=weeks[(endpoint[2] + 1):endpoint[3]],density=30,angle=0,col="dodgerblue2")
lines(x=weeks[(endpoint[1] + 1):endpoint[3]],cal_event_quantile1[(endpoint[1] + 1):endpoint[3]],lty = 3,lwd = 3, col="moccasin")
lines(x=weeks[(endpoint[2] + 1):endpoint[3]],cal_event_quantile3[(endpoint[2] + 1):endpoint[3]],lty = 4,lwd = 3,col="dodgerblue2")
legend("topleft",legend = c("Forecast 95% CI after 74 weeks","Forecast 95% CI after 77 weeks","Forecast 0.85 Quantile after 74 weeks","Forecast 0.85 Quantile after 77 weeks"), fill = c("moccasin", "dodgerblue2", NA, NA), density = c(30,30,0,0), angle = c(0,90,NA,NA), border = c("moccasin", "dodgerblue2", NA, NA), lty = c(NA,NA,3,4), lwd = c(NA, NA, 3,3), col = c(NA,NA, "moccasin", "dodgerblue2"),cex = 1, bty = "o", bg = "white" )
dev.off()


# We now turn to figure . While the forecast intervals after 74 and 78 weeks both accommodate the surge in deaths associated with the delta variant, the actual course of deaths was regarded as being much more plausible after 78 weeks than after 74 weeks, despite the fact that deaths had declined in that four week period. Evidently, the model picked up on the stalling disippation of deaths as evidence that transmission rate was actually increasing (which is also demonstrated in the concentration of the delta variant's coefficient). This is demonstrated even more strikingly in the change in forecast distribution for California
prior_rgb <- col2rgb('moccasin') / 255
post_rgb <- col2rgb('dodgerblue2') / 255
post2_rgb <- col2rgb('salmon4') / 255

png("output/us_analysis/endpoints/delta2_current_posterior.png", width = 6, height = 6, units = 'in', res = 300)
delta_effect_prior <- rnorm(10000, mean = 0, sd = sqrt(sfit[[names(endpoint[2])]]$prior$V0[2]))
hist(delta_effect_prior, col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75), breaks = seq(-20,20,0.2), freq = F, main = "Delta Prevalence Effect", ylim = c(0,3), xlim = c(-5,5), xlab = "Current Posterior", density = 20, angle=0)
hist(sfit[[names(endpoint[2])]]$samples$Xi0[,2], col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75), breaks = seq(-5,5,0.2), freq = F, add = T, density = 20, angle=45)
hist(sfit[[names(endpoint[3])]]$samples$Xi0[,2], col = rgb(post2_rgb[1], post2_rgb[2], post2_rgb[3], 0.75), breaks = seq(-5,5,0.2), freq = F, add = T, density = 20, angle=90)
legend("topleft", c("Prior - N(0,10)", "End of week 77", "End of week 90"), fill = c(rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 0.75),rgb(post_rgb[1], post_rgb[2], post_rgb[3], 0.75),rgb(post2_rgb[1], post2_rgb[2], post2_rgb[3], 0.75)), cex = 1, density = c(20,20,20), angle=c(0,45,90))
dev.off()

# Delta Prevalence plots
usa_delta_prevalence <- delta_prevalence %*% state_populations/sum(state_populations)
png("output/us_analysis/endpoints/usa_delta_prevalence.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:endpoint[3]], usa_delta_prevalence[PLOTSTART:endpoint[3]], xaxt = 'n', xlab = NA, ylab = "Proportion of New Cases", main = "USA Delta Variant Prevalence", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)
abline(v = weeks[endpoint[1]], col="moccasin")
abline(v = weeks[endpoint[2]], col = "dodgerblue2")
abline(v = weeks[endpoint[3]], col = "salmon4")
points(x = weeks[endpoint[1]], y = usa_delta_prevalence[endpoint[1]], pch = 4)
text(x = weeks[endpoint[1]], y = usa_delta_prevalence[endpoint[1]], labels = weeks[endpoint[1]], pos = 2)
points(x = weeks[endpoint[2]], y = usa_delta_prevalence[endpoint[2]], pch = 4)
text(x = weeks[endpoint[2]], y = usa_delta_prevalence[endpoint[2]], labels = weeks[endpoint[2]], pos = 4)
dev.off()

png("output/us_analysis/endpoints/cal_delta_prevalence.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:endpoint[3]], delta_prevalence[PLOTSTART:endpoint[3],k], xaxt = 'n', xlab = NA, ylab = "Proportion of New Cases", main = "CA Delta Variant Prevalence", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)
abline(v = weeks[endpoint[1]], col="moccasin")
abline(v = weeks[endpoint[2]], col = "dodgerblue2")
abline(v = weeks[endpoint[3]], col = "salmon4")
points(x = weeks[endpoint[1]], y = delta_prevalence[endpoint[1],k], pch = 4)
text(x = weeks[endpoint[1]], y = delta_prevalence[endpoint[1],k], labels = weeks[endpoint[1]], pos = 2)
points(x = weeks[endpoint[2]], y = delta_prevalence[endpoint[2],k], pch = 4)
text(x = weeks[endpoint[2]], y = delta_prevalence[endpoint[2],k], labels = weeks[endpoint[2]], pos = 4)
dev.off()

