library(smesir)
library(lubridate)
library(DescTools)
set.seed(123)

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



## --- FORECAST COMPARISON ---
PLOTSTART <- 50 # lower bound on time
prior_rgb <- col2rgb('moccasin') / 255
post_rgb <- col2rgb('dodgerblue2') / 255

# First version - parallel quantiles
qlevels <- seq(0.2,0.8,0.2)
event_quantiles74 <- apply(apply(sfor[[names(endpoint[1])]]$event_samples, c(1,2), sum),
                           MARGIN = 1,
                           FUN = quantile,
                           prob = qlevels)
event_quantiles77 <- apply(apply(sfor[[names(endpoint[2])]]$event_samples, c(1,2), sum),
                           MARGIN = 1,
                           FUN = quantile,
                           prob = qlevels)

png("output/us_analysis/endpoints/usa_endpoint_quantiles_forecasts.png", width = 8, height = 8, units = 'in', res = 300)
plot(weeks[PLOTSTART:endpoint[3]], rowSums(covid_deaths[PLOTSTART:endpoint[3],]), 
     ylim=c(0,2.6e4), xlim=c(weeks[PLOTSTART], weeks[endpoint[3]] + weeks(1)), xaxt = 'n', xlab = NA, ylab = NA, main = "US Weekly COVID-19 Deaths", type = "l",lwd=2)
axis(1, monthly,format(monthly, "%b %y"), cex.axis = 1)


for(r in 1:length(qlevels)){
  lines(x=weeks[(endpoint[1] + 1):endpoint[3]], 
        y = event_quantiles74[r,(endpoint[1] + 1):endpoint[3]],
        col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 1),
        lty = 3, lwd = 4)
  text(x=weeks[endpoint[3]],
       y = event_quantiles74[r,endpoint[3]],
       col = rgb(prior_rgb[1], prior_rgb[2], prior_rgb[3], 1),
       labels = as.character(qlevels[r]),
       adj = -0.1,
       cex = 1.2)
  lines(x=weeks[(endpoint[2] + 1):endpoint[3]], 
        y = event_quantiles77[r,(endpoint[2] + 1):endpoint[3]],
        col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 1),
        lty = 4, lwd = 4)
  text(x=weeks[endpoint[3]],
       y = event_quantiles77[r,endpoint[3]],
       col = rgb(post_rgb[1], post_rgb[2], post_rgb[3], 1),
       labels = as.character(qlevels[r]),
       adj = -0.1,
       cex = 1.2)
}
legend("topleft",legend = c("Forecast Quantiles after 74 weeks","Forecast Quantiles after 77 weeks"), lty = c(3,4), lwd = c(4,4), col = c("moccasin", "dodgerblue2"), cex = 1, bty = "o", bg = "white" )
dev.off()

## --- COVARIATE EFFECT POSTERIOR CONCENTRATION ---
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


## --- LOCAL COVARIATE EFFECTS AT 90 WEEKS ---
local_effect_quantiles <- apply(sfit[[names(endpoint[3])]]$samples$Xi[,2,], 2, quantile, probs=c(0.025, 0.5, 0.975))
local_effect_quantiles <- cbind(local_effect_quantiles, global = quantile(sfit[[names(endpoint[3])]]$samples$Xi0[,2], probs = c(0.025, 0.5, 0.975)))

png("output/us_analysis/endpoints/usa_local_effects.png", width = 6, height = 12, units = 'in', res = 300)
PlotDot(x = local_effect_quantiles[2,], 
        args.errbars = list(from = local_effect_quantiles[1,], to = local_effect_quantiles[3,]),
        xlim = c(0,2),
        main = "Local Delta Variant Effects")
dev.off()
